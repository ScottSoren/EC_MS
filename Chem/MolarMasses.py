import re

Zdict = {'H':1,'He':2,
	'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
	'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
	'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,
	'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,
	'Cs':55,'Ba':56,'La':57,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,
	'Fr':87,'Ra':88,'Ac':89,
	'Ce':58,'Gd':64,
	'U':92,
	'Heavy':0,'Light':0, 'e':0, '+':0, '-':0}
	# Heavy and Light are so that I can have isotopically labled stuff.
    # + and - are so that charges are collected by break_down

Mdict = {'H':1.008,'He':'4.002',
	'Li':6.941,'Be':9.012,'B':10.811,'C':12.011,'N':14.007,'O':15.999,'F':18.998,'Ne':20.180,
	'Na':22.990,'Mg':24.305,'Al':26.982,'Si':28.086,'P':30.974,'S':32.065,'Cl':35.453,'Ar':39.948,
	'K':39.098,'Ca':40.078,'Sc':44.956,'Ti':47.867,'V':50.942,'Cr':51.996,'Mn':54.938,'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.380,'Ga':69.723,'Ge':72.640,'As':74.922,'Se':78.96,'Br':79.904,'Kr':83.904,
	'Rb':85.468,'Sr':87.620,'Y':88.906,'Zr':91.224,'Nb':92.906,'Mo':95.960,'Tc':98,'Ru':101.070,'Rh':102.906,'Pd':106.420,'Ag':107.868,'Cd':112.411,'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.600,'I':126.904,'Xe':131.293,
	'Cs':132.905,'Ba':137.327,'La':138.905,'Hf':178.490,'Ta':180.948,'W':183.840,'Re':186.207,'Os':190.230,'Ir':192.217,'Pt':195.084,'Au':196.967,'Hg':200.590,'Tl':204.383,'Pb':207.200,'Bi':208.980,'Po':209,'At':210,'Rn':222,
	'Fr':223,'Ra':226,'Ac':227,
	'Ce':140.116,'Gd':157.250,
	'U':238.029,
	'Heavy':1,'Light':-1, # So that I can have isotopically labled stuff.
    'e':0, '-':0, '+':0}         # So that charges are collected by break_down

def BreakDown(*args, **kwargs):			#Reads parentheses and breaks down an arbitrary compound into a dictionary, like AxByCz to {'A':x,'B':y,'C':z}
    print('function Chem.BreakDown is now called Chem.break_down. Remember that next time!')
    return break_down(*args, **kwargs)

def break_down(compound, forgiving=True):
    '''
    Breaks a string representing a chemical formula down into constituent parts.
    Things in parentheses are considered one part.
    Any string of a capital letter followed by lower case letters is considered to be
        an irriducible element. 
    Any number is considered to quantify the element imediately proceeding it.
    Space is ignored
    Other characters raise a ValueError unless forgiving=True

    Example:
        >>> break_down('CH3(CH2)5CHO') 
        >>> {'C':2, 'H':4', 'CH2':5, 'O':1}
    
    This function is called recursively by get_elements.
    '''
    parts = {}
    number = ''
    element = ''
    subcompound = ''
    nest = 0
    N = len(compound)
    addit = False
    for i, char in enumerate(compound):
        if char is '(':
            if nest == 0:
                addit = True
            else:
                subcompound += char
            nest+=1
        elif nest > 0:
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
            elif re.search('[A-Z\-\+]',char):
                addit = True			
            elif re.search('\S',char):
                print('Not quite sure what you\'re talking about, mate, when you say ',char)
                if not forgiving:
                    raise ValueError        
        if i == N - 1:
            addit = True
        #print('char = ' + char + '\nelement = ' + element + '\naddit = ' + str(addit))       
        if addit:
            
            if len(number)>0:
                try:
                    n = int(number)
                except ValueError:
                    n = float(number)
                number = ''
            else:
                n = 1
            if len(element) > 0:
                if element in parts:
                    parts[element] += n
                else:
                    parts[element] = n
            if nest == 0:
                element = char
                if i == N - 1 and re.search('[A-Z\-\+]', char):
                    if element in parts:
                        parts[element] += 1
                    else:
                        parts[element] = 1
            addit = False
		
    return parts
	

def Mass(*args, **kwargs):
    print('Function Chem.Mass() has been replaced by Chem.get_mass(). Remember that next time!')
    return get_mass(*args, **kwargs)

def get_mass(compound, forgiving=True):			#Returns the molar mass of any chemical formula.

	if compound in Mdict:
		return(Mdict[compound])
	
	parts = break_down(compound, forgiving)
	
	M = 0
	for element in parts:
	#	A = input('checking for ' + element + ' in ' + Compound)
		if element==compound:
			print('Dude, man,', compound, 'just isn\'t a thing!')
			if not forgiving:
				raise ValueError
			return 0
		M += parts[element]*get_mass(element, forgiving=forgiving)
		
	return M

	
def SimpleForm(*args, **kwargs):
    print('function Chem.SimpleForm() is now called Chem.get_elements(). Remember that next time!!!')
    return get_elements(*args, **kwargs)

def get_elements(formula, forgiving=True):			
    #Determines number of each element in a complex formula.
    '''
    Given a chemical formula like 'CH3(CH2)5CHO', returns the number of each 
    element, here {'C':7, 'H':14, 'O':1}
    Works recursively if there are parentheses in the formula
    If forgiving=True, considers unidentifiable things to be elements
    '''
    elements = {} # the broken down formula will go here with elements as keys

    if type (formula) is str:  # this function breaks down each key of a dict. 
        formula = {formula:1}  # so for a string, just treat it as the one key of a dict
        
    cleared = []  # used to avoid an infinite loop in the case that something can't be broken down
    
    for compound in formula:
        nComp = formula[compound]

        if compound in elements:  #if it's already in elements, just add to the count!
            elements[compound] += nComp
            continue
        elif compound in Zdict: #if it's an irriducible atom, just put it in elements!
            elements[compound] = nComp 
            continue
    		#otherwise, we've got to break it down:
        parts = break_down(compound, forgiving=forgiving)
        
        if compound in parts and compound not in Zdict:  
            #If we can't break it down and it's not a recognizeable atom:
            print('dude,', compound, 'just isn\'t a thing')
            if not forgiving:
                raise ValueError
            cleared += [compound] # So that it doesn't put the function in an infinite loop
		
        for part, nPart in parts.items():  # Just everything from the breakdown into elements
            if part in ['g','l','s','aq']:
                continue # so that the state doesn't count as an element for e.g., CH3OH(aq)
            if part in elements:
                elements[part] += nPart * nComp
            else:
                elements[part] = nPart * nComp

	# Check if we put anything reducible into elements, in which case we've got to call this again.
    more = len([comp for comp in elements.keys() if comp not in Zdict and comp not in cleared])
    if more: 
	#	print('Simplifying', Output)
        return get_elements(elements)
	
    return elements



