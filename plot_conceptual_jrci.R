library(ggplot2)

ggplot() +
  # Random null model
  geom_line(data = data.frame(RR_target = seq(1, 100, 1)),
            aes(x = RR_target, y = 1 / 100),
            inherit.aes = FALSE,
            linewidth = 1,
            color = "blue") +
  annotate("text", x = 1, y = 0.05,
           hjust = "left",
           label = "Random null model",
           size = 6,
           color = "darkblue") + 
  # Identical systems model
  geom_line(data = data.frame(RR_target = seq(1, 100, 0.1)),
            aes(x = RR_target, y = 1 / RR_target),
            inherit.aes = FALSE,
            linewidth = 1,
            color = "blue") +
  annotate("text", x = 1.1, y = 1,
           hjust = "left",
           label = "Identical systems model",
           size = 6,
           color = "darkblue") + 
  geom_line(data = data.frame(RR_target = seq(1, 100, 1)),
            aes(x = RR_target, y = 0),
            inherit.aes = FALSE,
            linetype = "dashed",
            linewidth = 1,
            color = "black") +
  scale_x_log10() +
  xlab("Sub-system RR (%)") +
  ylab(TeX("$JRR / RR^2$")) +
  theme_classic() +
  theme(text = element_text(size = 22))

ggsave(filename = "plots/jrci_models.pdf",
       width = 1.5 * 4, height = 1.5 * 3)
