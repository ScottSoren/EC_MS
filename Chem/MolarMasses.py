
Zdict = {'H':1,'He':2,
	'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
	'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
	'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,
	'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,
	'Cs':55,'Ba':56,'La':57,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Ir':81,'Pb':82,'Bi':83,'Tl':84,'At':85,'Rn':86,
	'Fr':87,'Ra':88,'Ac':89,
	'Ce':58,'Gd':64,
	'U':92,
	'Heavy':0,'Light':0}
	#Heavy and Light are so that I can have isotopically labled stuff.

Mdict = {'H':1.008,'He':'4.002',
	'Li':6.941,'Be':9.012,'B':10.811,'C':12.011,'N':14.007,'O':15.999,'F':18.998,'Ne':20.180,
	'Na':22.990,'Mg':24.305,'Al':26.982,'Si':28.086,'P':30.974,'S':32.065,'Cl':35.453,'Ar':39.948,
	'K':39.098,'Ca':40.078,'Sc':44.956,'Ti':47.867,'V':50.942,'Cr':51.996,'Mn':54.938,'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.380,'Ga':69.723,'Ge':72.640,'As':74.922,'Se':78.96,'Br':79.904,'Kr':83.904,
	'Rb':85.468,'Sr':87.620,'Y':88.906,'Zr':91.224,'Nb':92.906,'Mo':95.960,'Tc':98,'Ru':101.070,'Rh':102.906,'Pd':106.420,'Ag':107.868,'Cd':112.411,'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.600,'I':126.904,'Xe':131.293,
	'Cs':132.905,'Ba':137.327,'La':138.905,'Hf':178.490,'Ta':180.948,'W':183.840,'Re':186.207,'Os':190.230,'Ir':192.217,'Pt':195.084,'Au':196.967,'Hg':200.590,'Tl':204.383,'Pb':207.200,'Bi':208.980,'Po':209,'At':210,'Rn':222,
	'Fr':223,'Ra':226,'Ac':227,
	'Ce':140.116,'Gd':157.250,
	'U':238.029,
	'Heavy':1,'Light':-1}
	#Heavy and light are so that I can have isotopically labled stuff.

def BreakDown(Compound, forgiving=True):			#Reads parentheses and breaks down an arbitrary compound into a dictionary, like AxByCz to {'A':x,'B':y,'C':z}
	import re
	Parts = {}
	number = ''
	element = ''
	subcompound = ''
	nest = 0
	N = len(Compound)
	i=1
	addit=0
	for char in Compound:
		if char is '(':
			if nest==0:
				addit=1
			else:
				subcompound += char
			nest+=1
		elif nest>0:
			if char is ')':
				nest -= 1
			if nest == 0:
				element = subcompound
				subcompound = ''
			else:
				subcompound += char	
		else:
			if re.search('[/.0-9]',char):
				number += char
			elif re.search('[a-z]',char):
				element += char	
			elif re.search('[A-Z]',char):
				addit = 1;			
			elif re.search('\S',char):
				print('Not quite sure what you\'re talking about, mate, when you say ',char)
				if not forgiving:
					raise ValueError        
		if addit == 1 or i == N:
			if len(number)>0:
				n = float(number)
				number = ''
			else:
				n=1
			if len(element)>0:
				if element in Parts:
					Parts[element]+=n
				else:
					Parts[element]=n
			if nest == 0:
				element = char
				if i==N and re.search('[A-Z]',char):
					if element in Parts:
						Parts[element]+=1
					else:
						Parts[element]=1
			addit=0
		i=i+1
		
	return Parts
	
def Mass(Compound, forgiving=True):			#Returns the molar mass of any chemical formula.

	if Compound in Mdict:
		return(Mdict[Compound])
	
	Parts = BreakDown(Compound, forgiving)
	
	M=0
	for element in Parts:
	#	A = input('checking for ' + element + ' in ' + Compound)
		if element==Compound:
			print('Dude, man,', Compound, 'just isn\'t a thing!')
			if not forgiving:
				raise ValueError
			return 0
		M += Parts[element]*Mass(element, forgiving=forgiving)
		
	return M
	
def SimpleForm(Input):			#Determines number of each element in a complex formula.
	
	cleared = ()
	from copy import deepcopy
		
	if type (Input) is str:
		Input = {Input:1}
	
	Output = deepcopy(Input)
	
	for Compound in Input:
		#print('In loop', Compound)
		if Compound in Zdict:
			continue
		#print('Still in loop')
		nComp = Input[Compound]
		Output[Compound]=0
		
		Parts = BreakDown(Compound)
		if Compound in Parts:
			print('dude,', Compound, 'just isn\'t a thing')
			cleared += (Compound,)
		
		for element in Parts:
			if element in Output:
				Output[element] += Parts[element]*nComp
			else:
				Output[element] = Parts[element]*nComp
				
		#print(Compound,'\n',Parts)
	
	Output = {i:Output[i] for i in Output if Output[i]!=0}
	
	more = 0
	for Compound in Output:
		#print('checking',Compound)
		if Compound in Zdict or Compound in cleared:
			continue
		else:
	#		print(Compound,'not okay')
			more +=1	
	if more>0: 
	#	print('Simplifying', Output)
		return SimpleForm(Output)
	
	return Output