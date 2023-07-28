from Ensemble.sparcc import SparCC, MakeBootstraps, PseudoPvals, EmpiricalBrownsMethod
import sys
import numpy

# Version of the Ensemble Plugin, but only runs SparCC
class SparCCPlugin:
   def input(self, filename):
      # Compute SparCC, Pearson and Spearman correlations on the input file (now just SparCC)
      #self.dirs = ["sparcc", "pearson", "spearman"]
      self.dirs = ["sparcc"]
      for algorithm in self.dirs:
         sys.argv = ["SparCC.py", filename, '-i', '5', "--cor_file=cor_"+algorithm+".out", "-a", algorithm]
         SparCC.driver()
      #   SparCC.driver(filename, ("-i", "5", "--cor_file=cor_"+algorithm+".out", "-a", algorithm))

      # Generate 100 permutations of the input file
      sys.argv = ["MakeBootstraps.py", filename, "-n", "1000", "-t", "permutation_#.txt"]
      MakeBootstraps.driver()

   def run(self):
      # Run SparCC, Pearson, and Spearman correlations on all permutations
      for algorithm in self.dirs:
         for i in range(1000):
            filename = "permutation_"+str(i)+".txt"
            sys.argv = ["SparCC.py", filename, "-i", "5", "--cor_file=perm_cor_"+algorithm+"_"+str(i)+".txt", "-a", algorithm]#, "--pval_file=pval_"+algorithm]
            SparCC.driver()

      # Compute unified p-values
      for algorithm in self.dirs:
         filename = "cor_"+algorithm+".out"
         filename2 = "perm_cor_"+algorithm+"_#.txt"
         sys.argv = ["PseudoPvals.py", filename, filename2, "1000", "-o" "pvals."+algorithm+".txt", "-t", "two_sided"]
         PseudoPvals.driver()


      using = 'sparcc'
      numsamples = 1

      # Quick read just to get number of OTUs
      tmpfilename = "perm_cor_sparcc_0.txt"
      tmpfile = open(tmpfilename, 'r')
      firstline = tmpfile.readline().strip()
      OTUs = firstline.split('\t')
      #OTUs.remove('OTU_id')
      self.numnodes = len(OTUs)
      tmpfile.close()

      # Initialize adjacency and sign matrices
      self.ADJ = []
      #SIGNS = []
      for i in range(self.numnodes):
           self.ADJ.append([])
           #SIGNS.append([])
           for j in range(self.numnodes):
              self.ADJ[i].append(0)#[j] = value
              #SIGNS[i].append("")
      
      # Adjacency matrix will be produced by SparCC (using)
      # However, any entries where all three algorithms
      # do not agree on sign will be set to zero
      for dir in self.dirs:
         myfile = "cor_" + dir + ".out"
         filestuff = open(myfile, 'r')
         firstline = filestuff.readline()
         lines = []
         for line in filestuff:
            lines.append(line)

         numlines = len(lines)
         self.bacteria = firstline.strip().split('\t')
         #self.bacteria = firstline.split('\t')
         #self.bacteria.remove('OTU_id')
         numbac = len(self.bacteria)

         for i in range(numlines):
            #print(lines[i])
            contents = lines[i].split('\t')
            if contents.count('\n') != 0:
               contents.remove('\n')
               contents.append('0.0')
            for j in range(numbac):
               value = float(contents[j+1].strip())
               if (dir == using):
                  self.ADJ[i][j] = value
               #if (value < 0):
               #   SIGNS[i][j] += "-"
               #elif (value > 0):
               #   SIGNS[i][j] += "+"
         #for i in range(self.numnodes):
         #    for j in range(self.numnodes):
         #       if self.ADJ[i][j] == 2:  # Now do one final sweep, setting flags to zero
         #         if (SIGNS[i][j].count('+') != 0 and SIGNS[i][j].count('-') != 0): #Disagree on sign
         #            self.ADJ[i][j] = 0
         filestuff.close()


         # Store correlations and pvalues for all three algorithms
         #correlations = dict()
         #pvalues = dict()
         #for dir in self.dirs:
         #    correlations[dir] = numpy.zeros([numsamples, self.numnodes, self.numnodes])
         #    pvalues[dir] = numpy.zeros([self.numnodes, self.numnodes])
         self.pvalues = numpy.zeros([self.numnodes, self.numnodes])


         #for dir in self.dirs:
         #samplefile = open("cor_" + dir + ".out")
         #sample = 0
         #samplefile.readline()
         #for i in range(self.numnodes):
         #line = samplefile.readline().strip()
         #line = line.split('\t')
         #line.remove(line[0])
         #for j in range(self.numnodes):
         #correlations[dir][sample][i][j] = line[j]
         dir = 'sparcc'
         pvaluefile = open("pvals."+str(dir)+".txt")
         pvaluefile.readline()
         for i in range(self.numnodes):
                pline = pvaluefile.readline().strip()
                pline = pline.split('\t')
                pline.remove(pline[0])
                for j in range(self.numnodes):
                      self.pvalues[i][j] = pline[j]



   def output(self, filename):
          filestuff2 = open(filename, 'w')
          filestuff2.write("\"\",")
          
          p_thresh = 0.01
          if (filename.endswith("unthresholded.csv")):
              p_thresh = 1
          for i in range(self.numnodes):
            for j in range(self.numnodes):
               if (i == j):
                   self.ADJ[i][j] = 1
                   self.ADJ[j][i] = 1
               elif (self.pvalues[i][j] > p_thresh):
                  #print(self.bacteria[i], self.bacteria[j], self.ADJ[i][j], self.pvalues[i][j])
                  self.ADJ[i][j] = 0
                  self.ADJ[j][i] = 0
               #else:
               #   print("PASSED:",self.bacteria[i],self.bacteria[j],self.ADJ[i][j], self.pvalues[i][j])


          for i in range(self.numnodes):
             filestuff2.write(self.bacteria[i])
             if (i != self.numnodes-1):
                filestuff2.write(",")
             else:
                filestuff2.write("\n")

          for i in range(self.numnodes):
                   filestuff2.write(self.bacteria[i]+',')
          #         filestuff2.write(self.samples[i]+',')
                   for j in range(self.numnodes):
                      filestuff2.write(str(self.ADJ[i][j]))
                      if (j < self.numnodes-1):
                         filestuff2.write(",")
                      else:
                         filestuff2.write("\n")

