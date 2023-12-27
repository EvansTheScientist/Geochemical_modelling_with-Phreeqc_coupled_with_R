

### Skript 11 # Projektarbeit T.Schnepper # 03.2021

## Custom Piper-Plot

# Credits: https://gist.github.com/johnDorian/5561272 # Jason Lessels jlessels@gmail.com


### This now consists of two functions. transform_piper_data transforms the data to match 
### the coordinates of the piper diagram. ggplot_piper does all of the background.


transform_piper_data <- function(Mg, Ca, Cl,SO4, name=NULL){
  if(is.null(name)){
    name = rep(1:length(Mg),3)
  } else {
    name = rep(name,3)
  }
  y1 <- Mg * 0.86603
  x1 <- 100*(1-(Ca/100) - (Mg/200))
  y2 <- SO4 * 0.86603
  x2 <-120+(100*Cl/100 + 0.5 * 100*SO4/100)
  new_point <- function(x1, x2, y1, y2, grad=1.73206){
    b1 <- y1-(grad*x1)
    b2 <- y2-(-grad*x2)
    M <- matrix(c(grad, -grad, -1,-1), ncol=2)
    intercepts <- as.matrix(c(b1,b2))
    t_mat <- -solve(M) %*% intercepts
    data.frame(x=t_mat[1,1], y=t_mat[2,1])
  }
  np_list <- lapply(1:length(x1), function(i) new_point(x1[i], x2[i], y1[i], y2[i]))
  npoints <- do.call("rbind",np_list)
  data.frame(observation=name,x=c(x1, x2, npoints$x), y=c(y=y1, y2, npoints$y))
}


ggplot_piper <- function() {
  library(ggplot2)
  grid1p1 <<- data.frame(x1 = c(20,40,60,80), x2= c(10,20,30,40),y1 = c(0,0,0,0), y2 = c(17.3206,34.6412,51.9618, 69.2824))
  grid1p2 <<- data.frame(x1 = c(20,40,60,80), x2= c(60,70,80,90),y1 = c(0,0,0,0), y2 = c(69.2824, 51.9618,34.6412,17.3206))
  grid1p3 <<- data.frame(x1 = c(10,20,30,40), x2= c(90,80,70,60),y1 = c(17.3206,34.6412,51.9618, 69.2824), y2 = c(17.3206,34.6412,51.9618, 69.2824))
  grid2p1 <<- grid1p1
  grid2p1$x1 <- grid2p1$x1+120
  grid2p1$x2 <- grid2p1$x2+120
  grid2p2 <<- grid1p2
  grid2p2$x1 <- grid2p2$x1+120
  grid2p2$x2 <- grid2p2$x2+120
  grid2p3 <<- grid1p3
  grid2p3$x1 <- grid2p3$x1+120
  grid2p3$x2 <- grid2p3$x2+120
  grid3p1 <<- data.frame(x1=c(100,90, 80, 70),y1=c(34.6412, 51.9618, 69.2824, 86.603), x2=c(150, 140, 130, 120), y2=c(121.2442,138.5648,155.8854,173.2060))
  grid3p2 <<- data.frame(x1=c(70, 80, 90, 100),y1=c(121.2442,138.5648,155.8854,173.2060), x2=c(120, 130, 140, 150), y2=c(34.6412, 51.9618, 69.2824, 86.603))
  
  p <- ggplot() +
    ## left hand ternary plot
    geom_segment(aes(x=0,y=0, xend=100, yend=0)) +
    geom_segment(aes(x=0,y=0, xend=50, yend=86.603)) +
    geom_segment(aes(x=50,y=86.603, xend=100, yend=0)) +
    ## right hand ternary plot
    geom_segment(aes(x=120,y=0, xend=220, yend=0)) +
    geom_segment(aes(x=120,y=0, xend=170, yend=86.603)) +
    geom_segment(aes(x=170,y=86.603, xend=220, yend=0)) +
    ## Upper diamond
    geom_segment(aes(x=110,y=190.5266, xend=60, yend=103.9236)) +
    geom_segment(aes(x=110,y=190.5266, xend=160, yend=103.9236)) +
    geom_segment(aes(x=110,y=17.3206, xend=160, yend=103.9236)) +
    geom_segment(aes(x=110,y=17.3206, xend=60, yend=103.9236)) +
    ## Add grid lines to the plots
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid1p1, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid1p2, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid1p3, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid2p1, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid2p2, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid2p3, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid3p1, linetype = "dashed", size = 0.5, colour = "grey50") +
    geom_segment(aes(x=x1, y=y1, yend=y2, xend=x2), data=grid3p2, linetype = "dashed", size = 0.5, colour = "grey50") +
    ### Labels and grid values
    #geom_text(aes(50,-10, label="Ca^2"), parse=T, size=4) + # Commented out, as parse=TRUE can cause issues
    
    geom_text(aes(c(20,40,60,80),c(-5,-5,-5,-5), label=c(80, 60, 40, 20)), size=5, face = "bold") +
    geom_text(aes(c(35,25,15,5),grid1p2$y2, label=c(80, 60, 40, 20)), size=5, face = "bold") +
    geom_text(aes(c(95,85,75,65),grid1p3$y2, label=c(80, 60, 40, 20)), size=5, face = "bold") +
    # geom_text(aes(17,50, label="Mg^2"), parse=T, angle=60, size=4) +
    coord_equal(ratio=1)+  
    geom_text(aes(17,50, label="Mg^'2+'"), angle=60, size=5, parse=TRUE, face = "bold") +  
    geom_text(aes(82.5,50, label="Na^+phantom()~+~K^+phantom()"), angle=-60, size=5, parse=TRUE, face = "bold") +
    geom_text(aes(50,-10, label="Ca^'2+'"), size=5, parse=TRUE, face = "bold") +
    
    
    geom_text(aes(170,-10, label="Cl^-phantom()"), size=5, parse=TRUE, face = "bold") +
    geom_text(aes(205,50, label="SO[4]^'2-'"), angle=-60, size=5, parse=TRUE, face = "bold") +
    geom_text(aes(137.5,50, label="CO[3]^'2-'~+~HCO[3]^-phantom()"), angle=60, size=5, parse=TRUE, face = "bold") +
    geom_text(aes(72.5,150, label="SO[4]^'2-'~+~Cl^-phantom()"), angle=60, size=5, parse=TRUE) +
    geom_text(aes(147.5,150, label="Ca^'2+'~+~Mg^'2+'"), angle=-60, size=5, parse=TRUE) + 
    
    geom_text(aes(c(155,145,135,125),grid2p2$y2, label=c(20, 40, 60, 80)), size=5, face = "bold") +
    geom_text(aes(c(215,205,195,185),grid2p3$y2, label=c(20, 40, 60, 80)), size=5, face = "bold") +
    geom_text(aes(c(140,160,180,200),c(-5,-5,-5,-5), label=c(20, 40, 60, 80)), size=5, face = "bold") +
    geom_text(aes(grid3p1$x1-5,grid3p1$y1, label=c(80, 60, 40, 20)), size=5, face = "bold") +
    geom_text(aes(grid3p1$x2+5,grid3p1$y2, label=c(20, 40, 60, 80)), size=5, face = "bold") +
    geom_text(aes(grid3p2$x1-5,grid3p2$y1, label=c(20, 40, 60, 80)), size=5, face = "bold") +
    geom_text(aes(grid3p2$x2+5,grid3p2$y2, label=c(80, 60, 40, 20)), size=5, face = "bold") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  return(p)
}

df_mgl <- Tidy_data_Greece


## Test-Daten simulieren

# Zwei Proben, Ionenkonzentrationen liegen mit Einheit mg/l vor

#df_mgl <- data.frame("Sample" = c("1", "2"),
#                     "Mg" = c(42.4, 27.0),
#                     "Ca" = c(92.0, 71.0),
#                     "Na" = c(45.5, 32.1),
#                     "K" = c(3.7, 6.4),
#                     "HCO3" = c(330.7, 330.7),
#                     "Cl" = c(138.65, 75.1),
#                     "SO4" = c(51.23, 64.76),
#                     "pH" = c(7.9, 8.2),
#                     "TDS" = c(688.8, 502.4))

#  meq/l aus mg/l berechnen: Massenkonzentration / Molare Masse * Valenz

df_meql <- df_mgl %>% summarize("Date" = Date,
                                "ID" = ID,  
                                "Mg_meq" = round(Mg * 2, digits = 2),
                                "Ca_meq" = round(Ca * 2, digits = 2),
                                "Na_meq" = round(Na * 1, digits = 2),
                                "K_meq" = round(K * 1, digits = 2),
                                "HCO3_meq" = round(HCO3 * 1, digits = 2),
                                "Cl_meq" = round(Cl * 1, digits = 2),
                                "SO4_meq" = round(SO4 * 2, digits = 2),
                                "pH" = pH) %>%
  filter(complete.cases(.) & ID == "ΥΝP_222" | 
           complete.cases(.) & ID == "MWS"| 
           complete.cases(.) & ID == "RW")

df_meql$ID_ext <- paste(df_meql$ID, df_meql$Date)

#  % meq/l aus meq/l berechnen

data <- df_meql %>% summarize("ID" = ID,
                              "ID_ext" = ID_ext,
                              "Ca" = round(Ca_meq / (Ca_meq + Mg_meq + Na_meq + K_meq) * 100, digits = 2),
                              "Mg" = round(Mg_meq / (Ca_meq + Mg_meq + Na_meq + K_meq) * 100, digits = 2),
                              "Cl" = round(Cl_meq/ (Cl_meq + SO4_meq + HCO3_meq) * 100, digits = 2),
                              "SO4" = round(SO4_meq / (Cl_meq + SO4_meq + HCO3_meq) * 100, digits = 2))

#transform the data into piper based coordinates
piper_data <- transform_piper_data(Ca=data$Ca, Mg = data$Mg, Cl=data$Cl, SO4= data$SO4, name=data$ID)
# The piper function now just plots the background
ggplot_piper()
# Now points can be added like...
ggplot_piper() + geom_point(aes(x,y), data=piper_data)
# colouring the points can be done using the observation value.
ggplot_piper() + geom_point(aes(x,y, colour= factor(observation)), data=piper_data)
# The size can be changed like..
ggplot_piper() + geom_point(aes(x,y, colour=factor(observation)), size=3, data=piper_data)
## Change colours and shapes and merging the legends together
ggplot_piper() + geom_point(aes(x,y, colour=factor(observation), fill=factor(observation), shape=factor(observation)), size=5, data=piper_data) + 
  scale_colour_manual(name="Samples", values=c("black","black","black"), labels=c("Mine water", "Rain water", "Groundwater")) +
  scale_fill_manual(name="Samples", values=c("#e8b56f", "#82d192", "#78aaee"), labels=c("Mine water", "Rain water", "Groundwater")) +
  scale_shape_manual(name="Samples", values=c(22,24,21), labels=c("Mine water", "Rain water", "Groundwater")) +
  theme(legend.position = c(0.9,0.55),
        legend.title = element_text(size = 15, family="verdana"),
        legend.text = element_text(size = 15, family="verdana"),
        axis.title.x = element_text(size = 15, family="verdana"),
        axis.title.y = element_text(size = 15, family="verdana"))
