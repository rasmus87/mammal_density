# Compare world mammal biomass numbers to Anthony D. Barnosky
# Sources: 
# https://www.pnas.org/doi/10.1073/pnas.0801918105
# https://ourworldindata.org/mammals

ggplot(mam.dens, aes(x = log10(10^log10BM/1000), y = log10.density.mean, col = Order.1.2, shape = Order.1.2)) +
  geom_point() +
  geom_smooth(data = mam.dens %>% filter(log10BM >= log10(44e3)),
              method = lm, aes(group = Order.1.2 == "Carnivora"), col = "black", lwd = .7) +
  annotate("rect", xmin = -4, xmax = log10(44), ymin = -3, ymax = 5,
           alpha = .5, fill = "grey") +
  theme_R() +
  geom_abline(slope = -0.44, intercept = 1.01, lty = 2) + ## Large herbivore density according to https://www.pnas.org/doi/10.1073/pnas.0801918105#sec-3
  geom_abline(slope = -1.31, intercept = 1.22, lty = 2) + ## Large carnivore density according to https://www.pnas.org/doi/10.1073/pnas.0801918105#sec-3
  geom_vline(xintercept = log10(44)) + ## Large carnivore density according to https://www.pnas.org/doi/10.1073/pnas.0801918105#sec-3
  xlab(expression(log[10]~Body~mass~(kg))) +
  ylab(expression(log[10]~Density~(km^-2))) +
  scale_color_manual(values = col27, breaks = orders, labels = labels) +
  scale_shape_manual(values = pch27, breaks = orders, labels = labels) +
  theme(legend.text = element_markdown()) +
  theme(legend.position = c(1, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(1,1)) +
  guides(col = guide_legend(nrow = 9)) +
  coord_cartesian(ylim = c(-1.5, 4.30), xlim = c(-2.7, 4))
ggsave("output/appendix1_fig1_BarnoskyComparison.png", width = 25.6, height = 16.5, units = "cm")

# Our results are very different than Barnosky (2008) since he does not use
# mammals below 44 kg, and use regressions like the dotted lines for the animal
# density which is very diffrent than our