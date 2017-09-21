library(easyPubMed)
library(tidyverse)
my_query <- '(1095-9203[TA] OR 0028-0836[TA] OR 0027-8424[TA]) AND 2016[DP]'
out <- batch_pubmed_download(pubmed_query_string = my_query, dest_file_prefix = "SNP", batch_size = 1000, dest_dir = '~/desktop/NHIR_textmining/')

out <- list.files(path = '~/desktop/NHIR_textmining/', pattern = "*.xml",full.names = T)

xmlTovec <- function(pubmed_data) {
  options(warn = -1)
  tmp.article <- custom_grep(xml_data = pubmed_data, tag = "PubmedArticle", format = "char")
  
  tmp.title <- custom_grep(xml_data = tmp.article, tag = "ArticleTitle", format = "char")
  tmp.abstract <- custom_grep(xml_data = tmp.article, tag = "AbstractText", format = "char")
  if (length(tmp.abstract) > 1){
    tmp.abstract <- paste(tmp.abstract, collapse = " ", sep = " ")
  } else if (length(tmp.abstract) < 1) {
    tmp.abstract <- NA
  }
  tmp.date <- custom_grep(xml_data = tmp.article, tag = "PubDate", format = "char")
  tmp.date <- sapply(c("Year", "Month", "Day"), (function(tt){
    custom_grep(xml_data = tmp.date, tag = tt, format = "char")   
  }))
  tmp.jabbrv  <- custom_grep(xml_data = tmp.article, tag = "ISOAbbreviation", format = "char")
  
  tmp.resout <- c(title=tmp.title,
                  abstract=tmp.abstract,
                  year = as.character(tmp.date[1]),
                  month = as.character(tmp.date[2]),
                  day = as.character(tmp.date[3]),
                  jabbrv=tmp.jabbrv)
  options(warn = 0)
  return(tmp.resout)
}
paperdf <- tibble()
for (i in 1:length(out)){
  xmlfile <- xmlParse(out[i])
  paper.data <- articles_to_list(xmlfile)
  papers.list <- t(sapply(1:length(paper.data), (function(i){
    out <-  tryCatch(xmlTovec(paper.data[i]), error = function(e) {NULL})
  })))
  papertemp <- as_tibble(papers.list)
  paperdf <- bind_rows(paperdf,papertemp)
}
paperdf

###

library(tidytext)
library(stringr)
library(ggplot2)

paperdf %>%
  group_by(jabbrv) %>%
  summarize(papers = n_distinct(title)) %>%
  ggplot(aes(jabbrv, papers)) +
  geom_col() +
  coord_flip()


##各期刊的前十大摘要高頻詞

wordf <- paperdf %>%
  mutate(line = 1:nrow(paperdf)) %>%
  filter(nchar(abstract) > 0) %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words) %>%
  filter(str_detect(word, "[^\\d]")) 

words_by_journal <- wordf %>%
  count(jabbrv, word, sort = TRUE) %>%
  ungroup()

tf_idf <- words_by_journal %>%
  bind_tf_idf(word, jabbrv, n) %>%
  arrange(desc(tf_idf))

tf_idf %>%
  group_by(jabbrv) %>%
  top_n(10, tf_idf) %>%
  ungroup() %>%
  mutate(word = reorder(word, tf_idf)) %>%
  ggplot(aes(word, tf_idf, fill = jabbrv)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ jabbrv, scales = "free") +
  ylab("tf-idf") +
  coord_flip()

##各期刊的前十大摘要高頻詞(以標題為主)

ordft <- paperdf %>%
  mutate(line = 1:nrow(paperdf)) %>%
  filter(nchar(title) > 0) %>%
  unnest_tokens(word, title) %>%
  anti_join(stop_words) %>%
  filter(str_detect(word, "[^\\d]")) 

words_by_journalt <- wordft %>%
  count(jabbrv, word, sort = TRUE) %>%
  ungroup()

tf_idft <- words_by_journalt %>%
  bind_tf_idf(word, jabbrv, n) %>%
  arrange(desc(tf_idf))

tf_idft %>%
  group_by(jabbrv) %>%
  top_n(10, tf_idf) %>%
  ungroup() %>%
  mutate(word = reorder(word, tf_idf)) %>%
  ggplot(aes(word, tf_idf, fill = jabbrv)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ jabbrv, scales = "free") +
  ylab("tf-idf") +
  coord_flip()

##詞關係

library(widyr)
library(igraph)
library(ggraph)

title_word_pairs <- wordft %>%
  pairwise_count(word,line,sort = TRUE)

set.seed(42)
title_word_pairs %>%
  filter(n >= 50) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Bigrams in title") +
  theme_void()

abs_word_pairs <- wordf %>%
  pairwise_count(word,line,sort = TRUE)

set.seed(42)
abs_word_pairs %>%
  filter(n >= 200) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Bigrams in abstract") +
  theme_void()

##

desc_dtm <- wordf %>%
  count(line, word, sort = TRUE) %>%
  ungroup() %>%
  cast_dtm(line, word, n)

library(topicmodels)
desc_lda <- LDA(desc_dtm, k = 20, control = list(seed = 42))
tidy_lda <- tidy(desc_lda)

top_terms <- tidy_lda %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)

top_terms %>%
  mutate(term = reorder(term, beta)) %>%
  group_by(topic, term) %>%    
  arrange(desc(beta)) %>%  
  ungroup() %>%
  mutate(term = factor(paste(term, topic, sep = "__"), 
                       levels = rev(paste(term, topic, sep = "__")))) %>%
  ggplot(aes(term, beta, fill = as.factor(topic))) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_x_discrete(labels = function(x) gsub("__.+$", "", x)) +
  labs(title = "Top 10 terms in each LDA topic",
       x = NULL, y = expression(beta)) +
  facet_wrap(~ topic, ncol = 5, scales = "free")

##

contributions <- wordf %>%
  inner_join(get_sentiments("afinn"), by = "word") %>%
  group_by(word) %>%
  summarize(occurences = n(),
            contribution = sum(score))

contributions %>%
  top_n(25, abs(contribution)) %>%
  mutate(word = reorder(word, contribution)) %>%
  ggplot(aes(word, contribution, fill = contribution > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip()

##
