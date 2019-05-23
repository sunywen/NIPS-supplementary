options(scipen=999)  # turn off scientific notation like 1e+06
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(quantmod)
aapl_data = getSymbols(Symbols = "AAPL", src = 'yahoo', auto.assign = FALSE, from = '2008-05-01', to = Sys.Date())
glimpse(aapl_data)

appl_data = read.csv('http://arn.la/applestock')

#ggplot(data=aapl_data, aes(x = AAPL.Open, y =  AAPL.Close)) + geom_line(color = "red")

ggplot(data = aapl_data, aes(x = date , y = AAPL.Close))+
  geom_line(color = "orange") + 
  scale_x_date(labels = date_format("%Y-%m"), breaks = date_breaks("6 months"))

g1 = g1 + theme_void()+ theme(axis.text = element_text(size=rel(0.8)))




