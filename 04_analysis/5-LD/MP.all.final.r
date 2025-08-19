pdf("test.all.final.pdf")
read.table("test.all.final.Japanese")->EJapanese;
plot(EJapanese[,1]/1000,EJapanese[,2],type="l",col="#00AFBB",main="LD decay",xlab="Distance(Kb)",xlim=c(0,300),ylim=c(0,0.521414665371257),ylab=expression(r^{2}),bty="n",lwd=2)
read.table("test.all.final.European")->EEuropean;
lines(EEuropean[,1]/1000,EEuropean[,2],col="#E7B800",lwd=2)

legend("topright",c("Japanese","European"),col=c("#00AFBB","#E7B800"),cex=1,lty=c(1,1),bty="n",lwd=2);
dev.off()
png("test.all.final.png",type="cairo")
read.table("test.all.final.Japanese")->EJapanese;
plot(EJapanese[,1]/1000,EJapanese[,2],type="l",col="red",main="LD decay",xlab="Distance(Kb)",xlim=c(0,300),ylim=c(0,0.521414665371257),ylab=expression(r^{2}),bty="n",lwd=2)
read.table("test.all.final.European")->EEuropean;
lines(EEuropean[,1]/1000,EEuropean[,2],col="black",lwd=2)

legend("topright",c("Japanese","European"),col=c("red","black"),cex=1,lty=c(1,1),bty="n",lwd=2);
dev.off()

