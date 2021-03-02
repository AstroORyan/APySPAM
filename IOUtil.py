"""
Edited: 02/02/2021.
Note, this script was edited by oryand on the above date to make the outputParticlesToFile function more general. The size_x and size_y functions
were updated with shapes. For original version, see the original JSPAM code.
"""
class IOUtil:

  @staticmethod
  def outputParticles(filename,x0):
    fo = open(filename,'w')
    IOUtil.outputParticlesToFile(fo,x0)


  @staticmethod
  def outputParticlesToFile(fo,x0):
    size_x = x0.shape[0]
    size_y = x0.shape[1]
    for i in range(size_x):
      dtmp = x0[i]
      for j in range(size_y):
        fo.write(IOUtil.formatDouble(dtmp[j]))
      # new line 
      fo.write("\n")
  
    fo.close()

  @staticmethod
  def formatDouble(num):
    return "%16.8e"%(num)
