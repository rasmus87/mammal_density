mam.dens.333 <- read_csv("builds/imputed.density_333.csv")
mam.dens.333$set <- "333"
mam.dens.33 <- read_csv("builds/imputed.density_33.csv")
mam.dens.33$set <- "33"
mam.dens.3 <- read_csv("builds/imputed.density_3.csv")
mam.dens.3$set <- "3"

mam.dens <- bind_rows(mam.dens.333, mam.dens.33, mam.dens.3)

str(mam.dens)

ggplot(mam.dens, aes(x = log10BM, y = log10.density.median, col = set)) +
  geom_point(pch = 1)

df <- data_frame(Binomial.1.2 = mam.dens.3$Binomial.1.2,
                 log10BM = mam.dens.3$log10BM,
                 Order.1.2 = mam.dens.3$Order.1.2,
                 med3 = mam.dens.3$log10.density.median,
                 med33 = mam.dens.33$log10.density.median,
                 med333 = mam.dens.333$log10.density.median,
                 sd3 = mam.dens.3$sd,
                 sd33 = mam.dens.33$sd,
                 sd333 = mam.dens.333$sd)

ggplot(df, aes(med3, med33)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1)
ggplot(df, aes(med3, med333)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1)
ggplot(df, aes(med33, med333)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1)

ggplot(df, aes(sd3, sd33)) +
  geom_point(pch = 1) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1)
ggplot(df, aes(sd3, sd333)) +
  geom_point(pch = 1) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1)
ggplot(df, aes(sd33, sd333)) +
  geom_point(pch = 1) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1)

df2 <- df %>% gather(key = "dataset", value = "sd", sd3:sd333)

ggplot(df2, aes(sd, col = dataset)) +
  geom_density()

