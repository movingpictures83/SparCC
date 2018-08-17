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
      sys.argv = ["MakeBootstraps.py", filename, "-n", "100", "-t", "permutation_#.txt"]
      MakeBootstraps.driver()

   def run(self):
      # Run SparCC, Pearson, and Spearman correlations on all permutations
      for algorithm in self.dirs:
         for i in range(100):
            filename = "permutation_"+str(i)+".txt"
            sys.argv = ["SparCC.py", filename, "-i", "5", "--cor_file=perm_cor_"+algorithm+"_"+str(i)+".txt", "-a", algorithm]
            SparCC.driver()

      # Compute unified p-values
      for algorithm in self.dirs:
         filename = "cor_"+algorithm+".out"
         filename2 = "perm_cor_"+algorithm+"_#.txt"
         sys.argv = ["PseudoPvals.py", filename, filename2, "100", "-o" "pvals."+algorithm+".txt", "-t", "two_sided"]
         PseudoPvals.driver()


      using = 'sparcc'
      numsamples = 1

      # Quick read just to get number of OTUs
      tmpfilename = "perm_cor_sparcc_0.txt"
      tmpfile = open(tmpfilename, 'r')
      firstline = tmpfile.readline().strip()
      OTUs = firstline.split('\t')
      OTUs.remove('OTU_id')
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
         self.bacteria = firstline.split('\t')
         self.bacteria.remove('OTU_id')
         numbac = len(self.bacteria)

         for i in range(numlines):
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
         pvalues = numpy.zeros([self.numnodes, self.numnodes])


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
                      pvalues[i][j] = pline[j]

         p_thresh = 0.01
         for i in range(self.numnodes):
            for j in range(self.numnodes):
               if (pvalues[i][j] > p_thresh):
                  self.ADJ[i][j] = 0
                  self.ADJ[j][i] = 0

         # Run Brown's Method, merge P-values
         #data_mat = numpy.zeros([len(self.dirs), 100])
         #p_values = numpy.zeros(len(self.dirs))


         # For each x and y:
         #   Read one element from each sample matrix
         #   Insert that element into the data matrix
         #   Get the p-values from the two-sided files, insert them into the pvalue vector
         #   Pass to the Brown function, the [0] element is the answer
         #merged_p_values = []

         #for x in range(self.numnodes):
         #   print "MERGING PVALUES FOR OTU ", OTUs[x]
         #   for y in range(self.numnodes):
         #      for dir in self.dirs:
         #         for sample in range(numsamples):
         #            data_mat[self.dirs.index(dir)][sample] = correlations[dir][sample][x][y]
         #         p_values[self.dirs.index(dir)] = pvalues[dir][x][y]
         #      if (y > x):
         #         merged_p_values.append((EmpiricalBrownsMethod.EmpiricalBrownsMethod(data_mat, p_values), (x, y)))

         # Hochberg Method for Cutting Off Pvalues
         #merged_p_values.sort()

         #m = len(merged_p_values)
         #q_star = 0.01   # What we use
         #for i in range(1, m+1):
         #   cutoff_P = (i/float(m))*q_star
         #   print "INDEX: ", i, " PVALUE: ", merged_p_values[i-1][0], " self.ADJACENCY VALUE: ", self.ADJ[merged_p_values[i-1][1][0]][merged_p_values[i-1][1][1]]
         #   if (merged_p_values[i-1][0] > cutoff_P):
         #      print "CUTOFF HIT AT INDEX ", i, " PVALUE: ", merged_p_values[i-1][0], " CUTOFF P: ", cutoff_P
         #      print "LENGTH OF MERGED VALUES: ", m
         #      for j in range(i+1, m+1):
         #         x = merged_p_values[j-1][1][0]
         #         y = merged_p_values[j-1][1][1]
         #         self.ADJ[x][y] = 0
         #         self.ADJ[y][x] = 0
         #      break


   def output(self, filename):
          filestuff2 = open(filename, 'w')
          filestuff2.write("\"\",")

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

