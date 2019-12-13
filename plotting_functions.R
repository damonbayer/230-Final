trace_plot <- function(data) {
  data %>% 
    pivot_longer(cols = everything()) %>% 
    group_by(name) %>% 
    mutate(t = row_number()) %>% 
    ggplot(aes(x = t, y = value, color = name)) +
    facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
    geom_line() +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle("Trace Plot")
  
}

dist_plot <- function(data) {
  data %>% 
    pivot_longer(cols = everything()) %>% 
    ggplot(aes(x = value, fill = name)) +
    facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
    geom_density() +
    theme_bw() +
    theme(legend.position="none") + 
    ggtitle("Posterior Distributions")  
}