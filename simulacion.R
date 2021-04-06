library(igraph)
library(rARPACK)
library(data.table)
library(magrittr)
library(ggplot2)


# -------------------------------------------------------------------------


algoritmo <- function(M, v2, k, c = 1.3, max_iter = 500){
  if (missing(v2)) {v2 <- eigs_sym(M, 2, "LA")$vectors[,2]}
  Q <- tail(order(abs(v2)),k)
  Q <- list(Q, sort(head(order(-rowSums(M[,Q])),k*c)))
  i <- 1
  while (length(intersect(Q[[i+1]],Q[[i]]))!=k && i < max_iter){
    i <- i+1
    # cat(i,' ')
    # flush.console()
    c2 <- rnorm(1,1,.3)
    if (k*c2 < 2) c2 <- c2 + 1 
    Q[[i+1]] <- sort(head(order(-rowSums(M[,Q[[i]]])), k*c2)) # aca hay problema si Q[[i]] es un solo elemento, por eso el -if- anterior
    # length(Q[[i]]) == 1
  }
  Q[[i+1]] <- sort(head(order(-rowSums(M[,Q[[i]]])),k))
  return(Q)
}



# -------------------------------------------------------------------------


n <- 1000
k <- 20
prob <- .95
P <- matrix(c(prob, .5, .5, .5), 2, 2)
grafo <- sample_sbm(n,P,c(k,n-k))

M <- as_adjacency_matrix(grafo,sparse=F)
v2 <- eigs_sym(M,2,"LA")$vectors[,2]


# -------------------------------------------------------------------------

set.seed(12345)


parametros <- setDT(expand.grid(k=20:50, n = 1000))
m_dentro <- 10
m_distinto <- 20


out <- lapply(split(parametros, by = c('n','k')), function(caso) {
  n <- caso$n
  k <- caso$k
  
  cat(n,k)
  replicate(m_distinto, {
  
  grafo <- sample_sbm(n,P,c(k,n-k))
  M <- as_adjacency_matrix(grafo,sparse=F)
  v2 <- eigs_sym(M,2,"LA")$vectors[,2]
  
  mean(replicate(m_dentro, {
    Q <- algoritmo(M, v2, k)
    identical(Q[[length(Q)]], 1:k)
  }))
  
  }, simplify = T)
    })

lapply(out, data.table) %>% 
  rbindlist(idcol = 'par') %>% 
  .[, c('n','k') := tstrsplit(par, '\\.')] %>% 
  .[,!'par'] %>% 
  .[, by = k, media := mean(V1)] %>%
  ggplot() + 
  geom_bin2d(aes(x = k, y = V1), drop = F) +
  geom_line(aes(x = k, y = media, group = 1, col = 'promedio')) + 
  ylab('P(Recuperar el subgrafo)') +
  scale_fill_gradient(low = "light blue", high = "dark blue") +
  scale_y_continuous(breaks = seq(0,1,.1)) +
  theme_bw() +
  ggtitle('n = 1000, p= 0.95') +
  ggsave('graf4.png')

