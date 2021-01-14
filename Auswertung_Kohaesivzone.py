# -*- coding: utf-8 -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy
from math import *
import mesh
import __main__
import regionToolset
import displayGroupMdbToolset as dgm
import xyPlot
import displayGroupOdbToolset as dgo
import random # Zufallszahlengenerator
import sys
import optimization
import os
import random
import shutil
import glob
##############################################################
class INPUTFILESEARCH:
  def __init__(self,inputfile):
    self.ifile		= inputfile
  ###
  def FindElsementSetForMaterial(self, question):
    mats	= []
    elsets	= []
    I 		= open(self.ifile,'r')
    endtrigger 	= 'ASSEMBLY'
    end 	= False
    line	= I.readline()
    while not end:
      if endtrigger in line:
	end 	= True
      if '*Solid Section' in line:
	mats.append(line.rstrip().split('material=')[1])
	elsets.append(line.split('elset=')[1].split(',')[0])
	line = I.readline()
      else:
	line = I.readline()
    I.close()
    #
    elsetname 	= elsets[mats.index(self.Choose(mats,question))]
    theset	= self.ReadElset(elsetname)
    theset = [ int(x) for x in theset ]
    return theset
  ###
  def ReadElset(self,name):
    I 		= open(self.ifile,'r')
    end 	= False
    line	= I.readline()
    while not end and line:
      if 'elset='+name in line:
	titleline 	= line
	end 		= True
	if 'generate' in titleline :
	  tmp		= I.readline()
	  startid	= tmp.rstrip().replace(' ','').split(',')[0]
	  endid		= tmp.rstrip().replace(' ','').split(',')[1]
	  idlist	= range(int(startid),int(endid)+1)
	else:
	  line 		= I.readline()
	  idlist	= []
	  while '*' not in line:
	    idlist+=line.rstrip().replace(' ','').split(',')
	    line = I.readline()
      else:
	line = I.readline()
    return idlist
  ###
  @staticmethod
  def Choose(choicelist,text):
    menu = []
    for c in choicelist:
      menu.append([c,'n'])
    selection = getInputs(menu, text)
    return choicelist[selection.index('j')]
  ###
  def CreateElementDict(self):
    I=open(self.ifile,'r')
    elementdict		= {}
    starttrigger 	= 'type=CPS4E'
    endtrigger		= 'set'
    start 		= False
    end 		= False
    line = I.readline()
    while not end and line:
      if start and endtrigger in line:
	end = True
      elif start and not end:
	element	= line.rstrip().replace(' ','').split(',')[0]
	n1	= line.rstrip().replace(' ','').split(',')[1]
	n2	= line.rstrip().replace(' ','').split(',')[2]
	n3	= line.rstrip().replace(' ','').split(',')[3]
	n4	= line.rstrip().replace(' ','').split(',')[4]
	elementdict[element] = [n1,n2,n3,n4]
	line = I.readline()
      elif starttrigger in line:
	start = True
	line = I.readline()
      else:
	line = I.readline()
    return elementdict
  ###
  def CreateNodeDict(self):
    I=open(self.ifile,'r')
    nodedict		= {}
    starttrigger 	= '*Node'
    endtrigger		= 'type=CPS4E'
    start 		= False
    end 		= False
    line = I.readline()
    while not end and line:
      if start and endtrigger in line:
	end = True
      elif start and not end:
	node	= line.rstrip().replace(' ','').split(',')[0]
	x	= line.rstrip().replace(' ','').split(',')[1]
	y	= line.rstrip().replace(' ','').split(',')[2]
	nodedict[node] = [x,y]
	line = I.readline()
      elif starttrigger in line:
	start = True
	line = I.readline()
      else:
	line = I.readline()
    return nodedict
##############################################################
class OutputCohesiveData:
  def __init__ (self, filename, cohesivedata, nodes, elements):
    print 'Schreibe Daten aus der Kohaesivzone'
    self.filename 	= filename
    self.data 		= cohesivedata
    self.nodes		= nodes
    self.elements	= elements
    self.upgradeelements()
    self.writedata()
   
  def upgradeelements(self):
    tmp = {}
    for e in self.elements.keys():
      x1 = float(self.nodes[self.elements[e][0]][0])
      x2 = float(self.nodes[self.elements[e][1]][0])
      x3 = float(self.nodes[self.elements[e][1]][0])
      x4 = float(self.nodes[self.elements[e][2]][0])
      tmp[str(e)] = {'nodes': self.elements[e],'nodecoordsX': [x1, x2, x3, x4], 'ipcoords': self.calcipcoords(self.elements[e])}
    self.elements = tmp

  def extrapolate(self, xGP1, xGP2, yGP1, yGP2, x1, x2):
    m = (yGP2 - yGP1)/(xGP2 - xGP1)
    b = yGP1 - m*xGP1
    return [m*x1+b, m*x2+b]
  
  def calcipcoords(self, elnodes):
    r1 = -1.0/sqrt(3)
    r2 =  1.0/sqrt(3)
    x = [float(self.nodes[str(elnodes[0])][0]),
	  float(self.nodes[str(elnodes[1])][0]),
	  float(self.nodes[str(elnodes[2])][0]),
	  float(self.nodes[str(elnodes[3])][0])]
    xIP1 = 0.5 * min(x) * (1.0 - r1) + 0.5 * max(x) * (1.0 + r1)
    xIP2 = 0.5 * min(x) * (1.0 - r2) + 0.5 * max(x) * (1.0 + r2)
    return [xIP1, xIP2]
  
  def writedata(self):
    ogp 	= open(self.filename+'_GP.gp' , 'w')
    onodes 	=  open(self.filename+'_Nodes.extra' , 'w')
    ogp.write('# Daten die an den Gausspunkten der Kohaesivelemente vorliegen\n')
    onodes.write('#Daten die von den Gausspunkten der Kohaesivelemente zu den Knoten extrapoliert wurden\n')
    ogp.write('# x	f-Damage	S12	S22	D_2	E_2\n')
    onodes.write('# x	f-Damage	S12	S22	D_2	E_2\n')
    xtmp= []
    for d in self.data.keys():
      xtmp.append(self.elements[d]['ipcoords'][0])
      
    
    for t in sorted(xtmp):
      for d in self.data.keys():  
	x1   =  min(self.elements[d]['nodecoordsX'])
	x2   =  max(self.elements[d]['nodecoordsX'])
	
	xIP1 =  self.elements[d]['ipcoords'][0]
	xIP2 =  self.elements[d]['ipcoords'][1]
	
	if xIP1 == t:
      
	  fIP1 = self.data[d]['Damage'][0]
	  fIP2 = self.data[d]['Damage'][1]
	  fnodes = self.extrapolate(xIP1, xIP2, self.data[d]['Damage'][0], self.data[d]['Damage'][1], x1, x2)
      
	  s12IP1 = self.data[d]['S12'][0]
	  s12IP2 = self.data[d]['S12'][1]
	  s12nodes = self.extrapolate(xIP1, xIP2, self.data[d]['S12'][0], self.data[d]['S12'][1], x1, x2)
      
	  s22IP1 = self.data[d]['S22'][0]
	  s22IP2 = self.data[d]['S22'][1]
	  s22nodes = self.extrapolate(xIP1, xIP2, self.data[d]['S22'][0], self.data[d]['S22'][1], x1, x2)
      
	  dIP1 = self.data[d]['D_2'][0]
	  dIP2 = self.data[d]['D_2'][1]
	  dnodes = self.extrapolate(xIP1, xIP2, self.data[d]['D_2'][0], self.data[d]['D_2'][1], x1, x2)
      
	  eIP1 = self.data[d]['D_2'][0]
	  eIP2 = self.data[d]['D_2'][1]
	  enodes = self.extrapolate(xIP1, xIP2, self.data[d]['E_2'][0], self.data[d]['E_2'][1], x1, x2)
      
	  ogp.write(str(xIP1)+ '	' + str(fIP1) + '	' + str(s12IP1) + '	' + str(s22IP1) + '	' + str(dIP1) + '	' + str(eIP1) + '\n')
	  ogp.write(str(xIP2)+ '	' + str(fIP2) + '	' + str(s12IP2) + '	' + str(s22IP2) + '	' + str(dIP2) + '	' + str(eIP2) + '\n')
	  ogp.write('#\n')
	  
	  onodes.write(str(x1)+ '	' + str(fnodes[0]) + '	' + str(s12nodes[0]) + '	' + str(s22nodes[0]) + ' 	' + str(dnodes[0]) + '	' + str(enodes[0]) + '\n') 
	  onodes.write(str(x2)+ '	' + str(fnodes[1]) + '	' + str(s12nodes[1]) + '	' + str(s22nodes[1]) + ' 	' + str(dnodes[1]) + '	' + str(enodes[1]) + '\n') 
    ogp.close()
    onodes.close()
##############################################################
def makeresultdir(name):
  filename = name.split('.')[0]+'_KoaesivEvaluation'
  if os.path.exists(filename):
    shutil.rmtree(filename)
    os.mkdir(filename)
  else:
    os.mkdir(filename)
  o = open('./'+filename+'/INFO.txt', 'w')
  o.write('Kurze Info zu den Dateiendungen: \n*.extra: Daten wurden von den Gausspunkten der Elemente zu den Knoten, aber nicht gemittelt. Dabei werden nur Daten aus der jeweiligen Section (lower/upper/cohesive ) berucksichtigt\n*.gp: Daten aus Gausspunkten\n*.nodal: Daten die direkt an den Knoten vorliegen (U,EPOT)')
  o.close()
  return filename
##############################################################
#theOdb 		= 'grosseDefo_UELS_v2.odb'
theOdb		= INPUTFILESEARCH.Choose(glob.glob('*odb'), 'Welcher Fall soll bearbeitet werden?')
outputdir = makeresultdir(theOdb)
odb = session.odbs[theOdb] # Festlegen der odb
theInstance = odb.rootAssembly.instances.keys()[-1] # Festlegen der Instance
AktInstance = odb.rootAssembly.instances[theInstance]
#
ifile 		= INPUTFILESEARCH(theOdb.replace('.odb','.inp'))
elements 	= ifile.CreateElementDict()
nodes 		= ifile.CreateNodeDict()
#
r 		= random.randint(1,10000)
cohesivesetname	= 'cohesive_elements'+str(r)
lowersetname	= 'lower_elements'+str(r)
uppersetname	= 'upper_elements'+str(r)
lowernsetname	= 'lower_nodes'+str(r)
uppernsetname	= 'upper_nodes'+str(r)

cohesive_eval = True
lower_el_eval = True
upper_el_eval = True
lower_nodal_eval = True
upper_nodal_eval = True

for __ in range(100):
  print ' '
################################################################
# Lese Daten aus Kohaesivelementen
################################################################
if cohesive_eval: # Nur zum Strukturieren
  cohesive_dict = {}
  cohesiveset	= ifile.FindElsementSetForMaterial('Welches Material ist in der Kohaesivzone?')
  for c in cohesiveset:
    cohesive_dict[str(c)] = {'Damage': [None, None, None, None], 'S12': [None, None, None, None], 'S22': [None, None, None, None], 'E_2': [None, None, None, None], 'D_2': [None, None, None, None]}
  fset = session.odbs[theOdb].rootAssembly.instances[theInstance].ElementSetFromElementLabels(name = cohesivesetname, elementLabels=tuple(cohesiveset))
  fField 	= odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['f-Damage'].getSubset(region = fset,position=INTEGRATION_POINT)
  for v in fField.values:
    cohesive_dict[str(v.elementLabel)]['Damage'][v.integrationPoint-1] = v.data
#
# Lese Spannungen aus den Kohaesivelementen
#  
  sField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['STRESS'].getSubset(region = fset,position=INTEGRATION_POINT)
  for v in sField.values:
    cohesive_dict[str(v.elementLabel)]['S22'][v.integrationPoint-1] = v.data[1]
    cohesive_dict[str(v.elementLabel)]['S12'][v.integrationPoint-1] = v.data[3]
#
# Lese D_2 aus den Kohaesivelementen
#  
  dField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['D-FIELD'].getSubset(region = fset,position=INTEGRATION_POINT)
  for v in dField.values:
    cohesive_dict[str(v.elementLabel)]['D_2'][v.integrationPoint-1] = v.data[1]
#
# Lese E_2 aus den Kohaesivelementen
#  
  eField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['E-FIELD'].getSubset(region = fset,position=INTEGRATION_POINT)
  for v in eField.values:
    cohesive_dict[str(v.elementLabel)]['E_2'][v.integrationPoint-1] = v.data[1]
  OutputCohesiveData('./'+outputdir+'/cohesive_elements',cohesive_dict, nodes, elements)
else:
  print 'Kohaesivdaten werden nicht geschrieben!'
################################################################
# Lese Daten aus den unteren Elementen
################################################################
def filterelements(elementset, nodeset, eldict):
  tmp = []
  nodeset = [int(n) for n in nodeset]
  for e in elementset:
    node1 = int(eldict[str(e)][0])
    node2 = str(eldict[str(e)][1])
    node3 = str(eldict[str(e)][2])
    node4 = str(eldict[str(e)][3])
    if node1 in nodeset:
      tmp.append(e)
    elif node2 in nodeset:
      tmp.append(e)
    elif node3 in nodeset:
      tmp.append(e)
    elif node4 in nodeset:
      tmp.append(e)
  return tmp
################################################################
y_coh = [float(nodes[str(elements[str(cohesiveset[0])][0])][1]),
    float(nodes[str(elements[str(cohesiveset[0])][1])][1]),
    float(nodes[str(elements[str(cohesiveset[0])][2])][1]),
    float(nodes[str(elements[str(cohesiveset[0])][3])][1])
    ]
y_lower = min(y_coh)
y_upper = max(y_coh)
tol = 1E-7
lower_nodes, upper_nodes = [],[]
for n in nodes.keys():
  if abs(float(nodes[n][1])-y_lower)<=tol:
    lower_nodes.append(int(n))
  elif abs(float(nodes[n][1])-y_upper)<=tol:
    upper_nodes.append(int(n))
# Suche die Knoten am unteren/oberen Ufer
if lower_el_eval: # Nur zum Strukturieren 
  lowernodedict = {}
  for l in lower_nodes:
    lowernodedict[str(l)] = {'x': nodes[str(l)][0], 'S12': [], 'S22': [], 'U11': [], 'U12': [], 'U22': [], 'U21': [], 'D_1': [], 'D_2': [], 'E_1': [], 'E_2': []}
  
  lowerelset	= ifile.FindElsementSetForMaterial('Welches Material ist unter der Kohaesivzone?')
  fset = session.odbs[theOdb].rootAssembly.instances[theInstance].ElementSetFromElementLabels(name = lowersetname, elementLabels=tuple(lowerelset))
  #
  sField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['STRESS'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in sField.values:
    if int(v.nodeLabel) in lower_nodes:
      lowernodedict[str(v.nodeLabel)]['S22'].append(v.data[1])
      lowernodedict[str(v.nodeLabel)]['S12'].append(v.data[3])
  # Uik11', 'Uik22', 'Uik33', 'Uik12'
  uField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['Uik'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in uField.values:
    if int(v.nodeLabel) in lower_nodes:
      lowernodedict[str(v.nodeLabel)]['U11'].append(v.data[0])
      lowernodedict[str(v.nodeLabel)]['U12'].append(v.data[3])
      lowernodedict[str(v.nodeLabel)]['U22'].append(v.data[1])
      lowernodedict[str(v.nodeLabel)]['U21'].append(v.data[2])

  dField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['D-FIELD'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in dField.values:
    if int(v.nodeLabel) in lower_nodes:
      lowernodedict[str(v.nodeLabel)]['D_1'].append(v.data[0])
      lowernodedict[str(v.nodeLabel)]['D_2'].append(v.data[1])
      
  eField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['E-FIELD'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in eField.values:
    if int(v.nodeLabel) in lower_nodes:
      lowernodedict[str(v.nodeLabel)]['E_1'].append(v.data[0])
      lowernodedict[str(v.nodeLabel)]['E_2'].append(v.data[1])
      
  s22 = open('./'+outputdir+'/S22_lower.extra', 'w')
  s12 = open('./'+outputdir+'/S12_lower.extra', 'w')
#            './'+outputdir+ /
  u11 = open('./'+outputdir+'/Uik11_lower.extra', 'w')
  u12 = open('./'+outputdir+'/Uik12_lower.extra', 'w')
  u21 = open('./'+outputdir+'/Uik21_lower.extra', 'w')
  u22 = open('./'+outputdir+'/Uik22_lower.extra', 'w')
#            './'+outputdir+ /
  d1  = open('./'+outputdir+'/D_1_lower.extra', 'w')
  d2  = open('./'+outputdir+'/D_2_lower.extra', 'w')
#            './'+outputdir+ /
  e1  = open('./'+outputdir+'/E_1_lower.extra', 'w')
  e2  = open('./'+outputdir+'/E_2_lower.extra', 'w')
# 
  for l in lowernodedict.keys():
    x = lowernodedict[l]['x']
    for s in lowernodedict[l]['S22']:
      s22.write(str(x) + '	' + str(s)+'\n')
    for s in lowernodedict[l]['S12']:
      s12.write(str(x) + '	' + str(s)+'\n') 
    for s in lowernodedict[l]['U11']:
      u11.write(str(x) + '	' + str(s)+'\n')
    for s in lowernodedict[l]['U12']:
      u12.write(str(x) + '	' + str(s)+'\n') 
    for s in lowernodedict[l]['U21']:
      u21.write(str(x) + '	' + str(s)+'\n') 
    for s in lowernodedict[l]['U22']:
      u22.write(str(x) + '	' + str(s)+'\n') 
    for s in lowernodedict[l]['D_1']:
      d1.write(str(x) + '	' + str(s)+'\n')
    for s in lowernodedict[l]['D_2']:
      d2.write(str(x) + '	' + str(s)+'\n') 
    for s in lowernodedict[l]['E_1']:
      e1.write(str(x) + '	' + str(s)+'\n') 
    for s in lowernodedict[l]['E_2']:
      e2.write(str(x) + '	' + str(s)+'\n') 
  s22.close()
  s12.close()
  u11.close()
  u12.close()
  u21.close()
  u22.close()
  d1.close()
  d2.close()
  e1.close()
  e2.close()
else:
  print 'Daten aus den unteren Kontinuumselementen werden nicht geschrieben!'
################################################################
# Lese Daten aus den oberen Elementen
################################################################
if upper_el_eval: # Nur zum Strukturieren 
  uppernodedict = {}
  for l in upper_nodes:
    uppernodedict[str(l)] = {'x': nodes[str(l)][0], 'S12': [], 'S22': [], 'U11': [], 'U12': [], 'U22': [], 'U21': [], 'D_1': [], 'D_2': [], 'E_1': [], 'E_2': []}
  upperelset	= ifile.FindElsementSetForMaterial('Welches Material ist ueber der Kohaesivzone?')
  fset = session.odbs[theOdb].rootAssembly.instances[theInstance].ElementSetFromElementLabels(name = uppersetname, elementLabels=tuple(upperelset))
  #
  sField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['STRESS'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in sField.values:
    if int(v.nodeLabel) in upper_nodes:
      uppernodedict[str(v.nodeLabel)]['S22'].append(v.data[1])
      uppernodedict[str(v.nodeLabel)]['S12'].append(v.data[3])
  # Uik11', 'Uik22', 'Uik33', 'Uik12'
  uField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['Uik'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in uField.values:
    if int(v.nodeLabel) in upper_nodes:
      uppernodedict[str(v.nodeLabel)]['U11'].append(v.data[0])
      uppernodedict[str(v.nodeLabel)]['U12'].append(v.data[3])
      uppernodedict[str(v.nodeLabel)]['U22'].append(v.data[1])
      uppernodedict[str(v.nodeLabel)]['U21'].append(v.data[2])

  dField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['D-FIELD'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in dField.values:
    if int(v.nodeLabel) in upper_nodes:
      uppernodedict[str(v.nodeLabel)]['D_1'].append(v.data[0])
      uppernodedict[str(v.nodeLabel)]['D_2'].append(v.data[1])
      
  eField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['E-FIELD'].getSubset(region = fset,position=ELEMENT_NODAL)
  for v in eField.values:
    if int(v.nodeLabel) in upper_nodes:
      uppernodedict[str(v.nodeLabel)]['E_1'].append(v.data[0])
      uppernodedict[str(v.nodeLabel)]['E_2'].append(v.data[1])
      
  s22 = open('./'+outputdir+'/S22_upper.extra', 'w')
  s12 = open('./'+outputdir+'/S12_upper.extra', 'w')
#            './'+outputdir+ /
  u11 = open('./'+outputdir+'/Uik11_upper.extra', 'w')
  u12 = open('./'+outputdir+'/Uik12_upper.extra', 'w')
  u21 = open('./'+outputdir+'/Uik21_upper.extra', 'w')
  u22 = open('./'+outputdir+'/Uik22_upper.extra', 'w')
#            './'+outputdir+ /
  d1  = open('./'+outputdir+'/D_1_upper.extra', 'w')
  d2  = open('./'+outputdir+'/D_2_upper.extra', 'w')
#            './'+outputdir+ /    upp
  e1  = open('./'+outputdir+'/E_1_upper.extra', 'w')
  e2  = open('./'+outputdir+'/E_2_upper.extra', 'w')
  
  for l in uppernodedict.keys():
    x = uppernodedict[l]['x']
    for s in uppernodedict[l]['S22']:
      s22.write(str(x) + '	' + str(s)+'\n')
    for s in uppernodedict[l]['S12']:
      s12.write(str(x) + '	' + str(s)+'\n') 
    for s in uppernodedict[l]['U11']:
      u11.write(str(x) + '	' + str(s)+'\n')
    for s in uppernodedict[l]['U12']:
      u12.write(str(x) + '	' + str(s)+'\n') 
    for s in uppernodedict[l]['U21']:
      u21.write(str(x) + '	' + str(s)+'\n') 
    for s in uppernodedict[l]['U22']:
      u22.write(str(x) + '	' + str(s)+'\n') 
    for s in uppernodedict[l]['D_1']:
      d1.write(str(x) + '	' + str(s)+'\n')
    for s in uppernodedict[l]['D_2']:
      d2.write(str(x) + '	' + str(s)+'\n') 
    for s in uppernodedict[l]['E_1']:
      e1.write(str(x) + '	' + str(s)+'\n') 
    for s in uppernodedict[l]['E_2']:
      e2.write(str(x) + '	' + str(s)+'\n') 
  s22.close()
  s12.close()
  u11.close()
  u12.close()
  u21.close()
  u22.close()
  d1.close()
  d2.close()
  e1.close()
  e2.close()
else:
  print 'Daten aus den oberen Kontinuumselementen werden nicht geschrieben!'
################################################################
# Lese Daten aus den unteren Knoten
################################################################
if lower_nodal_eval:
  lowerfacenodedict = {}
  for l in lower_nodes:
    lowerfacenodedict[str(l)] = {'x': None, 'U': [None, None], 'EPOT': None}
  
  lowerset = session.odbs[theOdb].rootAssembly.instances[theInstance].NodeSetFromNodeLabels(name = lowernsetname, nodeLabels=tuple(lower_nodes))
  
  xtmp = []
  uField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['U'].getSubset(region = lowerset, position=NODAL)
  for v in uField.values:
    xtmp.append(float(nodes[str(v.nodeLabel)][0]))
    lowerfacenodedict[str(v.nodeLabel)]['x'] = nodes[str(v.nodeLabel)][0]
    lowerfacenodedict[str(v.nodeLabel)]['U'][0] = v.data[0]
    lowerfacenodedict[str(v.nodeLabel)]['U'][1] = v.data[1]
  
  epotField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['EPOT'].getSubset(region = lowerset, position=NODAL)
  for v in epotField.values:
    lowerfacenodedict[str(v.nodeLabel)]['EPOT'] = v.data
  
  
  u = open('./'+outputdir+'/U_lower.nodal', 'w')
  e = open('./'+outputdir+'/EPOT_lower.nodal', 'w')
  u.write('#x	U_1	U_2')
  u.write('#x	EPOT')
  xtmp = sorted(xtmp)
  for xtest in xtmp:
    for l in lowerfacenodedict.keys():
      x = float(lowerfacenodedict[str(l)]['x'])
      if x == xtest:
	u.write(str(x)+'	'+str(float(lowerfacenodedict[str(l)]['U'][0]))+'	'+str(float(lowerfacenodedict[str(l)]['U'][1]))+'\n')
	e.write(str(x)+'	'+str(float(lowerfacenodedict[str(l)]['EPOT']))+'\n')
  u.close()
  e.close()
################################################################
# Lese Daten aus den oberen Knoten
################################################################
if upper_nodal_eval:
  upperfacenodedict = {}
  for l in upper_nodes:
    upperfacenodedict[str(l)] = {'x': None, 'U': [None, None], 'EPOT': None}
  
  upperset = session.odbs[theOdb].rootAssembly.instances[theInstance].NodeSetFromNodeLabels(name = uppernsetname, nodeLabels=tuple(upper_nodes))
  
  xtmp = []
  uField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['U'].getSubset(region = upperset, position=NODAL)
  for v in uField.values:
    xtmp.append(float(nodes[str(v.nodeLabel)][0]))
    upperfacenodedict[str(v.nodeLabel)]['x'] = nodes[str(v.nodeLabel)][0]
    upperfacenodedict[str(v.nodeLabel)]['U'][0] = v.data[0]
    upperfacenodedict[str(v.nodeLabel)]['U'][1] = v.data[1]
  
  epotField = odb.steps[odb.steps.keys()[-1]].frames[-1].fieldOutputs['EPOT'].getSubset(region = upperset, position=NODAL)
  for v in epotField.values:
    upperfacenodedict[str(v.nodeLabel)]['EPOT'] = v.data
  
  
  u = open('./'+outputdir+'/U_upper.nodal', 'w')
  e = open('./'+outputdir+'/EPOT_upper.nodal', 'w')
  u.write('#x	U_1	U_2')
  u.write('#x	EPOT')
  xtmp = sorted(xtmp)
  for xtest in xtmp:
    for l in upperfacenodedict.keys():
      x = float(upperfacenodedict[str(l)]['x'])
      if x == xtest:
	u.write(str(x)+'	'+str(float(upperfacenodedict[str(l)]['U'][0]))+'	'+str(float(upperfacenodedict[str(l)]['U'][1]))+'\n')
	e.write(str(x)+'	'+str(float(upperfacenodedict[str(l)]['EPOT']))+'\n')
  u.close()
  e.close()    