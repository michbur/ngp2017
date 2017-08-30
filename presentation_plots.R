library(ggplot2)
library(dplyr)

preds <- read.csv("negative_preds.csv")
pred_freq <- (as.vector(table(preds[["Probability"]] > 0.5)/nrow(preds))*100) %>% 
  formatC(digits = 2) %>% 
  paste0("%") %>% 
  data.frame(label = ., y = c(100, 50), x = c(0.25, 0.7), class = c("a", "b"))



cairo_ps(filename = "negative_hist.eps", width = 3*0.9, height = 3*0.9, pointsize = 24)
ggplot(preds, aes(x = Probability)) +
  geom_histogram(fill = rgb(0, 156, 219, maxColorValue = 255)) +
  geom_histogram(data = filter(preds, Probability > 0.5), 
                 fill = rgb(255, 65, 50, maxColorValue = 255)) +
  scale_x_continuous("Probability of amyloidogenicity") +
  scale_y_continuous("Number of non-amyloid peptides") +
  geom_text(data = pred_freq, aes(x = x, y = y, label = label, color = class), size = 10) +
  theme_bw() +
  scale_color_manual(values = c(rgb(0, 156, 219, maxColorValue = 255),
                                rgb(255, 65, 50, maxColorValue = 255))) +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(),
        plot.background = element_blank())
dev.off()



# peptides of interest --------------------------------
arrow_df <- data.frame(x = c(0.82, 0.82), y = c(70, 15))

cairo_ps(filename = "negative_hist_poi.eps", width = 3*0.9, height = 3*0.9, pointsize = 24)
ggplot(preds, aes(x = Probability)) +
  geom_histogram(fill = rgb(0, 156, 219, maxColorValue = 255)) +
  geom_histogram(data = filter(preds, Probability > 0.5), 
                 fill = rgb(255, 65, 50, maxColorValue = 255)) +
  scale_x_continuous("Probability of amyloidogenicity") +
  scale_y_continuous("Number of non-amyloid peptides") +
  geom_rect(xmin = 0.75, xmax = 0.9, ymin = -0.2, ymax = 10, color = "black", fill = NA) +
  geom_line(data = arrow_df, aes(x = x, y = y), arrow=arrow(length=unit(4, "mm"))) +
  theme_bw() +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(),
        plot.background = element_blank())
dev.off()

# results plot ------------------------------
ftir_res <- read.csv("ftir_amylogram_pasta2_literature.csv")

ftir_res_agg <- ftir_res %>% 
  group_by(amyl_db, amyl_ftir) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(amyl_db = factor(amyl_db, labels = c("Non-amyloid", "Amyloid")),
         amyl_ftir = factor(amyl_ftir, labels = c("Non-amyloid", "Amyloid"))) %>% 
  filter(amyl_db == "Non-amyloid",
         !is.na(amyl_ftir)) %>% 
  mutate(count = ifelse(count == 7, 4, count)) %>% 
  rbind(data.frame(amyl_db = "Non-amyloid", amyl_ftir = "Amyloid", count = 3)) %>% 
  mutate(frac = count/sum(count),
         count_nice = ifelse(count > 2, paste0(count, " peptides"), paste0(count, " peptide")),
         amyl_db = factor(c("Non-amyloid", 
                            "Amyloid", 
                            "Amyloid (literature)")),
         amyl_db = factor(amyl_db, levels = levels(amyl_db)[c(3, 1, 2)]))

cairo_ps(filename = "presentation_barlot.eps", width = 4, height = 2.7, pointsize = 24)
ggplot(ftir_res_agg, aes(x = amyl_ftir, y = frac, fill = amyl_db,
                         label = count_nice)) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5), color = "black", size = 4) +
  scale_y_continuous("Fraction of peptides") +
  scale_x_discrete("Experimental result") +
  scale_fill_manual("", values = c(rgb(0, 156, 219, maxColorValue = 255),
                                   rgb(255, 65, 50, maxColorValue = 255),
                                   rgb(255, 170, 0, maxColorValue = 255))) +
  theme_bw() +
  theme(legend.key.size = unit(0.3, "in"),
        legend.position = c(0.25, 0.75),
        legend.background = element_blank(),
        plot.background = element_blank())
dev.off()
