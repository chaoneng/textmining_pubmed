#Send search
library(RISmed)

res1 <- EUtilsSummary("data + science, population + health", 
                      type = "esearch", 
                      db = "pubmed",
                      datetype = "pdat",
                      retmax = 12000,
                      mindate = 2005, 
                      maxdate = 2016)

fetch <- EUtilsGet(res1, type = "efetch", db = "pubmed")

#Create data frame
abstracts <- data.frame(title = fetch@ArticleTitle,
                        abstract = fetch@AbstractText, 
                        journal = fetch@Title,
                        DOI = fetch@PMID, 
                        year = fetch@YearPubmed)
## ensure abstracts are character fields (not factors)
abstracts <- abstracts %>% mutate(abstract = as.character(abstract))
abstracts %>%
  head()

abstracts %>%
  group_by(year) %>%
  count() %>%
  filter(year > 2013) %>%
  ggplot(aes(year, n)) +
  geom_point() +
  geom_line() +
  labs(title = "Pubmed articles with search terms `data science` & `population health` \n2015-2016", hjust = 0.5,
       y = "Articles")

#Word cloud
cloud <- abstracts %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words) %>%
  count(word, sort = TRUE) 

cloud %>%
  with(wordcloud(word, n, min.freq = 10, max.words = 1000, colors = brewer.pal(8, "Dark2")), scale = c(8,.3), per.rot = 0.4)

#Bigram wordcloud
bigrams_united %>%
  with(wordcloud(bigram, n, max.words = 1000, random.order = FALSE, colors = brewer.pal(9, "Set1"), scale = c(8, 0.3)), per.rot = 0.4)

#Journal wordlcoud
cloud3 <- abstracts %>%
  select(journal) %>%
  group_by(journal) %>%
  count(sort = TRUE)
cloud3 %>%
  with(wordcloud(journal, n, min.freq = 10, random.order = FALSE, max.words = 80, colors = brewer.pal(9, "Set1")), rot.per = .6)

#Extract abstracts containing the phrase ‘data science’
g <- abstracts[grepl("data science", abstracts$abstract),]
g1 <- g$DOI %>% list
abstracts <- abstracts %>% 
  mutate(DOI  = as.character(DOI)) 

abstracts[abstracts$DOI %in% g1[[1]],] %>%
  select(title, journal, DOI) %>%
  knitr::kable()

#Extract abstracts containing the phrase ‘big data’
g <- abstracts[grepl("big data", abstracts$abstract),]
g1 <- g$DOI %>% list
abstracts <- abstracts %>% 
  mutate(DOI  = as.character(DOI)) 

abstracts[abstracts$DOI %in% g1[[1]],] %>%
  select(title, journal, DOI, year) %>%
  knitr::kable()

#Create document term matrix
abstracts %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words) %>%
  count(DOI, word, sort = TRUE) %>%
  cast_dtm(DOI, word, n) ->
  abstracts1

#Topic modelling (cluster analysis of abstracts)
library(topicmodels)
abs_lda <- LDA(abstracts1, k = 10, control = list(seed = 1234))
abs_lda_td <- tidytext:::tidy.LDA(abs_lda)

#Classify documents
abs_lda_gamma <- tidytext:::tidy.LDA(abs_lda, matrix = "gamma")
abs_class <- abs_lda_gamma %>%
  group_by(document) %>%
  top_n(1, gamma) %>%
  ungroup() %>%
  arrange(gamma)
abs_class %>%
  sample_n(6)

abstracts %>%
  group_by(journal) %>%
  count(sort = TRUE) %>%
  filter(n >=40) ->top40
abstracts %>%
  left_join(top40) %>%
  filter(!is.na(n)) %>%
  rename(document = DOI) %>%
  left_join(abs_class) %>%
  filter(!is.na(topic)) %>%
  ggplot(aes(document, factor(topic))) +
  geom_jitter(aes(colour = factor(topic)), size = 1, width = 0.1) +
  facet_wrap(~journal) +
  theme(strip.text.x = element_text(size = 8), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size  =10)) +
  scale_color_viridis(discrete = TRUE, option = "C") +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  labs(title = "Top 15 journals",
       subtitle = "Documents by topic",
       x = "Time")

