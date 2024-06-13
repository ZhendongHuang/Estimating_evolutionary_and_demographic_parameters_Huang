epsilon = 0.001
sam_size = 40

name = "gihatcon_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"
file = open(name,"rb")
gihatcon = pickle.load(file)
file.close()

name = "gitrue_" + "eps=" + str(epsilon) + "_sam_size=" + str(sam_size) + "constne_len108.dat"
file = open(name,"rb")
gitrue = pickle.load(file)
file.close()


k1= scipy.stats.gaussian_kde(gitrue, bw_method=None, weights=None)
k2= scipy.stats.gaussian_kde(gihatcon, bw_method=None, weights=None)

#ktg.set_bandwidth(bw_method=kerden.factor / 3.)

plt.hist(gihatcon, color = 'blue', edgecolor = 'black',bins=50,density=True,alpha=0.5)

xx= numpy.linspace(0,200000,1000)
y1 = k1.pdf(xx)
y2 = k2.pdf(xx)



plt.plot(xx,y1,'k-.')
plt.plot(xx,y2,'r')
plt.xlim(0,100000)
plt.ylim(0,5*10**(-5))

plt.legend(['True density','Estimated density','Histogram of the estimates'],fontsize=14,loc = 'upper right')


plt.xlabel('time to MRCA (generations)',fontsize=16)
plt.ylabel('density',fontsize=16)
plt.title('Time to MRCA of Model C',fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)