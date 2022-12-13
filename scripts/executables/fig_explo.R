source("scripts/functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("simul/param_explo", recursive = FALSE)
modulo <- pi
#####################

#Data
df.fig_explo <- df.data(sims.dirs, pattern = "simul/param_explo/", variable="netw", file_size=100000, w_of_6=TRUE, network=FALSE)
df.fig_explo$param <- str_split(df.fig_explo$data.dir, "/", simplify=TRUE)[,3]
df.fig_explo$value <-  str_split(df.fig_explo$param, "\\_",simplify=TRUE)[,2]
df.fig_explo$param <-  str_split(df.fig_explo$param, "\\_",simplify=TRUE)[,1]


################################################################################

values_lab <- c(
  `04` = "1",
  `06` = "2*",
  `12` = "8",
  `20` = "16",
  `30` = "26")

pfig <- ggplot(data=subset(df.fig_explo, param=="altgenes"), aes(corrS, corrM))+
  geom_point(aes(col=ecc_M), alpha=0.3, show.legend = TRUE)+scale_color_viridis_c(option = "plasma", limits=c(min(df.fig_explo$ecc_M), max(df.fig_explo$ecc_M)))+
  coord_fixed(ratio = 1, xlim = c(-0.9,0.9), ylim = c(-0.9,0.9), expand = TRUE, clip = "on")+
  labs(y="r(M)", x="r(S)", title="TF number", col = "e(M)")+
  facet_wrap(value ~., labeller = as_labeller(values_lab), ncol=5) + theme_bw()+
  theme(plot.margin = margin(t=4,0,0,0, "lines"),legend.direction="horizontal", legend.position = c(0.5, 1.5))
pfig

values_lab <- as_labeller(c(
  `00050` = "50",
  `01000` = "1000*",
  `00500` = "500",
  `05000` = "5000",
  `10000` = "10000"))

pfig1 <- ggplot(data=subset(df.fig_explo, param=="pop"), aes(corrS, corrM))+
  geom_point(aes(col=ecc_M), alpha=0.3, show.legend = FALSE)+
  labs(y="r(M)", x="r(S)", title="Population size" )+scale_color_viridis_c(option = "plasma", limits=c(min(df.fig_explo$ecc_M), max(df.fig_explo$ecc_M)))+
  coord_fixed(ratio = 1, xlim = c(-0.9,0.9), ylim = c(-0.9,0.9), expand = TRUE, clip = "on")+theme_bw()+
  facet_wrap(~value, labeller = values_lab, ncol=6)+
  theme(plot.margin = margin(t=0,0,0,0, "lines"))  #, strip.text = element_blank()
pfig1

values_lab <- as_labeller(c(
  `0.01` = "0.01",
  `0.1` = "0.1*",
  `0.5` = "0.5",
  `0.05` = "0.05",
  `0.9` = "0.9"))

pfig2 <- ggplot(data=subset(df.fig_explo, param=="mutsd"), aes(corrS, corrM))+
  geom_point(aes(col=ecc_M), alpha=0.3, show.legend = FALSE)+
  labs(y="r(M)", x="r(S)", title="Mutsd" )+scale_color_viridis_c(option = "plasma", limits=c(min(df.fig_explo$ecc_M), max(df.fig_explo$ecc_M)))+
  coord_fixed(ratio = 1, xlim = c(-0.9,0.9), ylim = c(-0.9,0.9), expand = TRUE, clip = "on")+
  facet_wrap(value ~., labeller = as_labeller(values_lab), ncol=6) + theme_bw() +
  theme(plot.margin = margin(t=0,0,0,0, "lines"))
pfig2

values_lab <- as_labeller(c(
  `0.1` = "0.1",
  `10` = "10*",
  `100` = "100",
  `1000` = "1000",
  `1` = "1"))

pfig3 <- ggplot(data=subset(df.fig_explo, param=="Ssize"), aes(corrS, corrM))+
  geom_point(aes(col=ecc_M), alpha=0.3, show.legend = FALSE)+
  labs(y="r(M)", x="r(S)", title="Selection strength" )+scale_color_viridis_c(option = "plasma", limits=c(min(df.fig_explo$ecc_M), max(df.fig_explo$ecc_M)))+
  coord_fixed(ratio = 1, xlim = c(-0.9,0.9), ylim = c(-0.9,0.9), expand = TRUE, clip = "on")+
  facet_wrap(value ~., labeller = as_labeller(values_lab), ncol=6) + theme_bw() +
  theme(plot.margin = margin(t=0,0,0,0, "lines"))
pfig3

values_lab <- as_labeller(c(
  `0.1` = "0.1",
  `0.5` = "0.5*",
  `0.25` = "0.25",
  `0.75` = "0.75",
  `0.9` = "0.9"))

pfig4 <- ggplot(data=subset(df.fig_explo, param=="basal"), aes(corrS, corrM))+
  geom_point(aes(col=ecc_M), alpha=0.3, show.legend = FALSE)+
  labs(y="r(M)", x="r(S)", title="Basal expression" )+scale_color_viridis_c(option = "plasma", limits=c(min(df.fig_explo$ecc_M), max(df.fig_explo$ecc_M)))+
  coord_fixed(ratio = 1, xlim = c(-0.9,0.9), ylim = c(-0.9,0.9), expand = TRUE, clip = "on")+
  facet_wrap(value ~., labeller = as_labeller(values_lab), ncol=6) + theme_bw() +
  theme(plot.margin = margin(t=0,0,0,0, "lines"))
pfig4

cairo_pdf("figures/fig_explo.pdf", width=9, height=12)
grid.arrange(
  pfig,pfig1,pfig2,pfig3,pfig4,
  ncol = 1,
  nrow = 5,
  clip = FALSE,
  heights = c(0.99,0.73,0.73,0.73,0.73)
)
dev.off()




################################################################################
sims.dirs <- c("simul/param_explo/mutsd_0.5", "simul/param_explo/mutsd_0.9")

df.G <- data.frame(NULL)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  df <- df.fig1(sims.dir, all.gen=TRUE)
  df.G <- rbind(df.G, df)
}
df.G[,9] <- str_split(df.G$data.dir, "/", simplify = TRUE)[,3]


ggplot(df.G[2748:2788,], aes(x = Gen, y=ang_M, col=V9))+geom_point()


