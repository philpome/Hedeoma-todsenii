library(ggmap)
library(maps)
library(mapdata)

#API key redacted for public repo: 
#ggmap::register_google(key="")

#map of NM region:
map <- get_googlemap(center = c(lon = -106.301024, lat = 33.114069),
                     zoom = 9, scale = 2,
                     maptype ='satellite',
                     color = 'color')
p<-ggmap(map)

#Polygons for mountain ranges:
SALongitude <- c(-106.8117687, -106.8200085, -106.7788097, -106.7458507, -106.6607067,
                 -106.6085216, -106.4986584, -106.3915417, -106.3723156, -106.4327404, 
                 -106.4904186, -106.4766857, -106.4162609, -106.3503429, -106.3778088,
                 -106.4931652, -106.6661999, -106.8117687)
SALatitude <- c(33.267406, 33.1731999, 32.9890916, 32.7999811, 32.5827001, 32.453005,
                32.441416, 32.503979, 32.6405392, 32.8138321, 32.9706595, 33.1364091,
                33.2421411, 33.3615105, 33.4463469, 33.4303031, 33.3982067, 33.267406)
SAdata <- as.data.frame(cbind(SALongitude,SALatitude))

SacLongitude <- c(-106.0070202, -106.0042736, -105.9932873, -105.949342, -105.9713146,
                  -105.949342, -105.9191296, -105.8229992, -105.6774303, -105.6224987,
                  -105.6609508, -105.6472179, -105.6884167, -105.7268688, -105.6884167,
                  -105.6966564, -105.699403, -105.8422253, -106.0070202)
SacLatitude <- c(33.4325952, 33.2582196, 33.1410088, 33.0213383, 32.9015052, 32.7653441,
                 32.6336005, 32.5364019, 32.552609, 32.654415, 32.8045983, 32.933784,
                 33.0512712, 33.1456083, 33.3018463, 33.4257186, 33.5219418, 33.5448363, 
                 33.4325952)
Sacdata <- as.data.frame(cbind(SacLongitude, SacLatitude))

#Integrate polygons into map:
library(ggspatial)
p2<-p+
  geom_polygon(data=SAdata,aes(x=SALongitude,y=SALatitude),alpha=0.3,colour="deepskyblue3",fill="deepskyblue3")+
  geom_polygon(data=Sacdata,aes(x=SacLongitude,y=SacLatitude),alpha=0.3,colour="coral1",fill="coral1")

#Calculate map min/max: 
x_lim <- p2$data[c(1, 4), 1] * c(1, .9998)
y_lim <- p2$data[c(1, 4), 2] * c(1.003, 1)
x_lim
y_lim

#Add scalebar: 
library(ggsn)
p3 <- p2 +
  scalebar(x.min = -107.0, x.max = -105.5,
           y.min = 32.4, y.max = 33.8,
           dist = 25, dist_unit="km", transform = TRUE, model = 'WGS84',
           location="topright", st.size=4, border.size=.5, st.color="white") +
  blank()
p3 
ggsave("hed-map.png",dpi=300)