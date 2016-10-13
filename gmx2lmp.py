import sys
top=sys.argv[1]       
xyz=sys.argv[2]
typelist=[]
if ".gro" in xyz:
 print  "gro"
 filetype = "gro"
 with open(xyz, 'r') as pos:
  xyzs = pos.readlines()
  atoms = int(xyzs[1])
  boundx,boundy,boundz = xyzs[atoms+2].split() 
  boundx=float(boundx)*10
  boundy=float(boundy)*10
  boundz=float(boundz)*10
  xyzs=xyzs[2:]
#  print xyzs[0]
elif ".pdb" in xyz:
 print "pdb"
 filetype = "pdb"
 with open(xyz, 'r') as pos:
  xyzs = pos.readlines()
  print xyzs[2]
  for start in range(0,20):
    if "CRYST" in xyzs[start]:
     boundx = xyzs[start][8:17]  
     boundy = xyzs[start][17:26]  
     boundz = xyzs[start][26:35]  
    elif "ATOM" in xyzs[start]:
     xyzs=xyzs[start:]
     break
    start=start+1
mols=-1
startatom=0
boundread = file("bound.tmp",'w')  
atomread = file("atom.tmp",'w')  
bondread = file("bond.tmp",'w')  
angleread = file("angle.tmp",'w')  
diheread = file("dihe.tmp",'w')  
coefread = file("coef.txt",'w')  
impread = file("imp.tmp",'w') 
atomlist=[]
masslist=[]
bondlist=[]
bondnums=0
anglelist=[] 
exanglelist=[] 
anglenums=0 
tmpangle=0 
dihelist=[] 
exdihelist=[] 
dihenums=0 
tmpdihe=0 
implist=[] 
eximplist=[] 
tmpimp=0 
impnums=0 
import subprocess
def grep(filename, arg):
    process = subprocess.Popen(['grep', '-n', arg, filename], stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    tmp= stdout.split(':',1)
    if  len(tmp)==1:
	line=0
    else:
        line,data= tmp
        line=int(line)
    return line

def getxyz(num):    # del1ed  atoms  set this to 99999 for no change
  if filetype == "gro":
    #a,b,c,x,y,z = xyzs[num].split()  
    x = xyzs[num][20:28]
    y = xyzs[num][28:36]  
    z = xyzs[num][36:44]  
    x=float(x)*10  
    y=float(y)*10  
    z=float(z)*10  
  elif filetype == "pdb":
    x = xyzs[num][30:38]  
    y = xyzs[num][38:46]  
    z = xyzs[num][46:54]  
  return  x,y,z 
x,y,z=getxyz(0)
with open(top, 'r') as fd:
    for line in fd:
        mol=line.rstrip( )
        molline=0
        mols += 1
        atomf=grep(mol,'atoms')
        bondf=grep(mol,'bonds') 
        pairf=grep(mol,'pairs') 
        anglef=grep(mol,'angles') 
        dihef=grep(mol,'dihedrals ] ; propers') 
        impf = grep(mol,'dihedrals ] ; impropers') 
        atomstart = atomf+2
        atomend = 0 
        bondstart = 0 
        bondend = 0 
        anglestart = 0 
        angleend = 0 
        dihestart = 0 
        diheend = 0 
        impstart = 0 
        impend = 0 
        with open(mol, 'r') as fe:
          endfile =len(fe.readlines()) 
        if bondf ==0:
            atomend = atomf+2 
        elif anglef ==0:
            atomend = bondf-1 
            bondstart = bondf+2  
            bondend = bondf+2 
        elif pairf ==0:
            atomend = bondf-1 
            bondstart = bondf+2  
            bondend = anglef-1 
            anglestart = anglef+2 
            angleend = endfile 
            dihestart = 0 
            diheend = 0 
            impstart = 0 
            impend = 0 
        elif impf ==0:
            atomend = bondf-1 
            bondstart = bondf+2  
            bondend = pairf-1 
            anglestart = anglef+2 
            angleend = dihef-1 
            dihestart = dihef+3  
            diheend = endfile 
            impstart = 0 
            impend = 0 
        else :
            atomend = bondf-1 
            bondstart = bondf+2  
            bondend = pairf-1 
            anglestart = anglef+2 
            angleend = dihef-1 
            dihestart = dihef+3  
            diheend = impf-1 
            impstart = impf+3 
            impend = endfile 
        with open(mol, 'r') as fe:
         for tmpline in fe:
          molline = molline + 1
          if(molline >= atomstart)&(molline <= atomend):
            sa = tmpline.split()
            nr,ntype,resnr,resid,atom,cgnr,charge,mass,t1,t2,t =sa   
            if ntype not in atomlist:
              atomlist.append(ntype)
              masslist.append(mass)
            nr=int(nr)
            x,y,z = getxyz(startatom+nr) 
            atomread.write("%7s  %7s  %7s %17s %17s %17s %17s\r\n"%(startatom+nr,mols,atomlist.index(ntype)+1,charge,x,y,z))
            if molline==atomend:
               nums=nr
          if(molline >= bondstart)&(molline <= bondend):
            sa = tmpline.split()
            at1,at2,funct,r,k,t,n1,t2,n2 =sa   
            if [funct,r,k] not in bondlist:
              bondlist.append([funct,r,k])
            at1=int(at1)
            at2=int(at2)
            bondnums=bondnums+1
            bondread.write("%7s  %7s  %7s %7s  #  %7s - %7s   \r\n"%(bondnums,bondlist.index([funct,r,k])+1,startatom+at1,startatom+at2,n1,n2))
          if(molline >= anglestart)&(molline <= angleend):
            sa = tmpline.split()
            at1,at2,at3,funct,r,k,t,n1,t2,n2,t3,n3 =sa   
            if [funct,r,k] not in anglelist:
              anglelist.append([funct,r,k])
            at1=int(at1)
            at2=int(at2)
            at3=int(at3) 
            anglenums=anglenums+1
            angleread.write("%7s  %7s  %7s  %7s  %7s  #  %7s - %7s - %7s  \r\n"%(anglenums,anglelist.index([funct,r,k])+1,startatom+at1,startatom+at2,startatom+at3,n1,n2,n3))
          if(molline >= dihestart)&(molline <= diheend):
            sa = tmpline.split()
            at1,at2,at3,at4,funct,r,k,o,t,n1,n2,n3,n4 =sa   
            if [funct,r,k,o] not in dihelist:
              dihelist.append([funct,r,k,o])
            at1=int(at1)
            at2=int(at2)
            at3=int(at3) 
            at4=int(at4) 
            dihenums=dihenums+1
            diheread.write("%7s  %7s  %7s  %7s  %7s %7s #  %7s - %7s - %7s -%7s \r\n"%(dihenums,dihelist.index([funct,r,k,o])+1,startatom+at1,startatom+at2,startatom+at3,startatom+at4,n1,n2,n3,n4))
            tmpdihe = [at1,at2,at3,at4]
          if(molline >= impstart)&(molline <= impend):
            sa = tmpline.split()
            at1,at2,at3,at4,funct,r,k,o,t,n1,n2,n3,n4 =sa   
            if [funct,r,k,o] not in implist:
                implist.append([funct,r,k,o])
            at1=int(at1)
            at2=int(at2)
            at3=int(at3) 
            at4=int(at4) 
            impnums=impnums+1
            impread.write("%7s  %7s  %7s  %7s  %7s %7s #  %7s - %7s - %7s -%7s \r\n"%(impnums,implist.index([funct,r,k,o])+1,startatom+at1,startatom+at2,startatom+at3,startatom+at4,n1,n2,n3,n4))
            tmpimp = [at1,at2,at3,at4]
          if(molline == endfile):
            startatom = startatom+ nums
          
boundread.write(" lammps data write from %17s & %17s \r\n"%(top,xyz))
boundread.write("%17s atoms \r\n"%(startatom ))
boundread.write("%17s bonds \r\n"%( bondnums ))
boundread.write("%17s angles \r\n"%( anglenums ))
boundread.write("%17s dihedrals \r\n "%( dihenums ))
boundread.write("%17s impropers \r\n"%( impnums ))
boundread.write("\r\n ")
boundread.write("%17s atom types \r\n"%( len(atomlist) ))
boundread.write("%17s bond types \r\n"%( len(bondlist) ))
boundread.write("%17s angle types \r\n"%( len(anglelist) ))
boundread.write("%17s dihedral types \r\n "%( len(dihelist) ))
boundread.write("%17s improper types \r\n"%( len(implist) ))
boundread.write("\r\n ")
boundread.write("%17f  %17f  xlo xhi \r\n"%(0,boundx ))
boundread.write("%17f  %17f  ylo yhi \r\n"%(0,boundy ))
boundread.write("%17f  %17f  zlo zhi \r\n"%(0,boundz ))
boundread.write("\r\n  Masses \r\n")
for ats in  atomlist:
     boundread.write("%7s %12s  # %7s  \r\n"%(atomlist.index(ats)+1, masslist[atomlist.index(ats)] ,ats ))
 
boundread.write("\r\n  Bonds \r\n\r\n")
for ats in  bondlist:
     boundread.write("%7s %12s  # bond  \r\n"%(bondlist.index(ats)+1, ats ))
boundread.write("\r\n Angles   \r\n\r\n")
for ats in  anglelist:
     boundread.write("%7s %12s  # angle \r\n"%(anglelist.index(ats)+1, ats ))
boundread.write("\r\n  Dihedrals \r\n\r\n")
for ats in  dihelist:
     boundread.write("%7s %12s  # dihe \r\n"%(dihelist.index(ats)+1, ats ))
boundread.write("\r\n  Impropers \r\n\r\n")
for ats in  implist:
     boundread.write("%7s %12s  # imp \r\n"%(implist.index(ats)+1, ats ))
boundread.write("\r\n ")




atomread.close()
boundread.close()
bondread.close()
angleread.close()
diheread.close()
