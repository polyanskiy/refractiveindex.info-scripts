###########################################################################################
#         THIS PROGRAM IS IN PUBLIC DOMAIN                                                #
#         COPYRIGHT AND RELATED RIGHTS WAIVED VIA CC0 1.0                                 #
###########################################################################################

#                 agf2yml (Zemax to refractiveindex.info converter)                       #

#-----------------------------------------------------------------------------------------#
#  author: Mikhail Polyanskiy (polyanskiy@refractiveindex.info)                           #
#  dependencies: Python3, NumPy                                                           #
#-----------------------------------------------------------------------------------------#

import os, math

#agffile = 'input/schottzemax-20150722.agf'
#ymldir = 'output/schott'
#references = '1) <a href=\\"http://refractiveindex.info/download/data/2015/schott-optical-glass-collection-datasheets-july-2015-us.pdf\\">SCHOTT optical glass data sheets 2015-07-22</a><br>2) <a href=\\"http://refractiveindex.info/download/data/2015/schottzemax-20150722.agf\\">SCHOTT Zemax catalog 2015-07-22</a>'

#agffile = 'input/OHARA_151201.agf'
#ymldir = 'output/ohara'
#references = '1) <a href=\\"http://refractiveindex.info/download/data/2015/ohara_2015-12-01.pdf\\">OHARA optical glass datasheets 2015-12-01</a><br>2) <a href=\\"http://refractiveindex.info/download/data/2015/OHARA_151201.agf\\">OHARA Zemax catalog 2015-12-01</a>'

#agffile = 'input/HIKARI.agf'
#ymldir = 'output/hikari'
#references = '1) <a href=\\"http://refractiveindex.info/download/data/2015/HIKARI_Catalog.pdf\\">HIKARI optical glass catalog 2015-04-01</a><br>2) <a href=\\"http://refractiveindex.info/download/data/2015/HIKARI.agf\\">HIKARI Zemax catalog</a>'

agffile = 'input/HOYA20150618.agf'
ymldir = 'output/hoya'
references = '<a href=\\"http://refractiveindex.info/download/data/2015/HOYA20150618.agf\\">HOYA Zemax catalog 2015-06-18</a>'

wl = []
IT = []


def ClearGlassData():
    name = ''
    comments = ''
    formula = 0
    glasscode = ''
    nd = 0
    vd = 0
    status = 0 # 1:Standard, 2:Preferred, 3:Special, 4:Obsolete, 5:Melt
    meltfreq = 0
 
    wlmin = 0
    wlmax = 0

    disp_formula_coefficients = ''

    thermal_disp_coefficients = ''

    cte1 = '-' #coefficient of termal expansion -30C - +70C [1E-6/K]
    cte2 = '-' #coefficient of termal expansion +20C - +300C [1E-6/K]
    density = '-' # [g/cm^3]
    dpgf = '-' #dPgF

    CR = '-' # Climatic resistance: 1(high) - 4(low)
    FR = '-' # Stain resistance: 0(high) - 5(low)
    SR = '-' # Acid resistance: 1(high) - 4(low) and 51-53(very low)
    AR = '-' # Alkali resistance: 1(high) - 4(low)
    PR = '-' # Phosphate resistance: 1(high) - 4(low)

    del wl[:]
    del IT[:]



def WriteYML(): 
    
    if glass_count == 1:
        if not os.path.exists(ymldir):
            os.mkdir(ymldir)
        os.chdir(ymldir)
        
    print(str(glass_count)+': '+name)
    
    try:
        ymlfile = open(name+'.yml', 'w+', encoding='utf-8')
    except TypeError: # Python2 has no encoding= argument
        ymlfile = open(name+'.yml', 'w+')

    ymlfile.write('# this file is part of refractiveindex.info database\n# refractiveindex.info database is in the public domain\n# copyright and related rights waived via CC0 1.0\n\n')
    ymlfile.write('REFERENCES: "' + references + '"\n')
    
    if comments:
        ymlfile.write('COMMENTS: "' + ' '.join(comments) + '"\n')
    
    ymlfile.write('DATA:\n')

    if formula == "1" or formula == "13":
        ymlfile.write('  - type: formula 3 \n')
    if formula == "2":
        ymlfile.write('  - type: formula 2 \n')
    ymlfile.write('    range: ' + format(float(wlmin)) + ' ' + format(float(wlmax)) + '\n')
    ymlfile.write('    coefficients:')
    if formula == "1" or formula == "13":
        if float(disp_formula_coefficients[0]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[0])))
        if float(disp_formula_coefficients[1]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[1])) + ' 2')
        if float(disp_formula_coefficients[2]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[2])) + ' 4')
        if float(disp_formula_coefficients[3]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[3])) + ' -2')
        if float(disp_formula_coefficients[4]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[4])) + ' -4')
        if float(disp_formula_coefficients[5]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[5])) + ' -6')
        if float(disp_formula_coefficients[6]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[6])) + ' -8')
        if float(disp_formula_coefficients[7]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[7])) + ' -10')
        if float(disp_formula_coefficients[8]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[8])) + ' -12')
    if formula == "2":
        ymlfile.write(' 0')
        if len(disp_formula_coefficients)>1 and float(disp_formula_coefficients[0]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[0])) + ' ' + format(float(disp_formula_coefficients[1])))
        if len(disp_formula_coefficients)>3 and  float(disp_formula_coefficients[2]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[2])) + ' ' + format(float(disp_formula_coefficients[3])))
        if len(disp_formula_coefficients)>5 and  float(disp_formula_coefficients[4]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[4])) + ' ' + format(float(disp_formula_coefficients[5])))
        if len(disp_formula_coefficients)>7 and  float(disp_formula_coefficients[6]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[6])) + ' ' + format(float(disp_formula_coefficients[7])))
        if len(disp_formula_coefficients)>9 and  float(disp_formula_coefficients[8]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[8])) + ' ' + format(float(disp_formula_coefficients[9])))
        if len(disp_formula_coefficients)>11 and  float(disp_formula_coefficients[10]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[10])) + ' ' + format(float(disp_formula_coefficients[11])))
        if len(disp_formula_coefficients)>13 and  float(disp_formula_coefficients[12]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[12])) + ' ' + format(float(disp_formula_coefficients[13])))
        if len(disp_formula_coefficients)>15 and  float(disp_formula_coefficients[14]):
            ymlfile.write(' ' + format(float(disp_formula_coefficients[14])) + ' ' + format(float(disp_formula_coefficients[15])))
    ymlfile.write('\n')
    
    ymlfile.write('  - type: tabulated k\n')
    ymlfile.write('    data: |\n')  
    for i in range (0,len(wl)):
        if float(IT[i]) != 0:
            k = -float(wl[i])/(4*math.pi)*math.log(float(IT[i]))/10000
            ymlfile.write('        ' + "%.3f"%(float(wl[i])) + ' ' + format(k, '.4E') + '\n')
    
    ymlfile.write('INFO:\n')
    ymlfile.write('    n_is_absolute: false\n')
    ymlfile.write('    λ_is_vacuum: false\n')
    if len(thermal_disp_coefficients)>6:
        ymlfile.write('    temperature: ' + format(float(thermal_disp_coefficients[6])) + ' °C\n')
    if len(thermal_disp_coefficients)>5 and (float(thermal_disp_coefficients[0]) or float(thermal_disp_coefficients[1]) or float(thermal_disp_coefficients[2]) or float(thermal_disp_coefficients[3]) or float(thermal_disp_coefficients[4]) or float(thermal_disp_coefficients[5])): 
        ymlfile.write('    coefficients_of_thermal_dispersion: ' + ' '.join([format(float(thermal_disp_coefficients[i])) for i in range(6) ]) + '\n') 
    ymlfile.write('    n_d: ' + nd + '\n')
    ymlfile.write('    V_d: ' + vd + '\n')
    if float(glasscode)>100000:
        ymlfile.write('    glass_code: ' + glasscode + '\n')
    if float(status)==1:    
        ymlfile.write('    glass_status: standard\n')
    if float(status)==2:    
        ymlfile.write('    glass_status: preferred\n')
    if float(status)==3:    
        ymlfile.write('    glass_status: special\n')
    if float(status)==4:    
        ymlfile.write('    glass_status: obsolete\n')
    if float(status)==5:    
        ymlfile.write('    glass_status: melt\n')
    if meltfreq != 0:
        ymlfile.write('    glass_melt_frequency: ' + meltfreq + '\n')
    if density != '-' and float(density):
        ymlfile.write('    density: ' + format(float(density)) + ' g/cm<sup>3</sup>\n')
    if cte1 != '-' and float(cte1): #-30...70 C
        ymlfile.write('    coefficient_of_thermal_expansion_1: ' + format(round(float(cte1)*1.0e-6,15)) + ' K<sup>-1</sup>\n')
    if cte2 != '-' and float(cte2): #20...300 C
        ymlfile.write('    coefficient_of_thermal_expansion_2: ' + format(round(float(cte2)*1.0e-6,15)) + ' K<sup>-1</sup>\n')
    if dpgf != '-' and float(dpgf):
        ymlfile.write('    ΔP_gF: ' + format(float(dpgf)) + '\n')
    if CR != '-' and float(CR) != -1:
        ymlfile.write('    climatic_resistance: ' + format(float(CR)) + '\n')
    if FR != '-' and float(FR) != -1:
        ymlfile.write('    stain_resistance: ' + format(float(FR)) + '\n')
    if SR != '-' and float(SR) != -1:
        ymlfile.write('    acid_resistance: ' + format(float(SR)) + '\n')
    if AR != '-' and float(AR) != -1:
        ymlfile.write('    alkali_resistance: ' + format(float(AR)) + '\n')
    if PR != '-' and float(PR) != -1:
        ymlfile.write('    phosphate_resistance: ' + format(float(PR)) + '\n')

    ymlfile.close()
    
    print('ok')

 


### main program
    
ClearGlassData()

glass_count = 0 #glass count

with open(agffile) as agf:
    for line in agf:
        data = line.split()
        if not len(data):
            continue
        if data[0] == 'NM':
            if glass_count!=0: WriteYML()
            ClearGlassData()
            glass_count +=1
            name = data[1]
            formula = data[2]
            glasscode = data[3]
            nd = data[4]
            vd = data[5]
            if len(data) >= 8: status = data[7] 
            if len(data) >= 9 and data.count('-') < 0: meltfreq = data[8]
            else: meltfreq = 0
            continue
        if data[0] == 'LD':
            wlmin = data[1]
            wlmax = data[2]
            continue
        if data[0] == 'CD':
            data.remove('CD')
            disp_formula_coefficients = data
            continue
        if data[0] == 'TD':
            data.remove('TD')
            thermal_disp_coefficients = data
            continue
        if data[0] == 'ED':
            cte1 = data[1]
            cte2 = data[2]
            density = data[3]
            dpgf = data[4]
            continue
        if data[0] == 'OD':
            #? = data[1]
            if len(data)>2:
                CR = data[2]
                if CR == '3-4':
                    CR = '3.5'
                FR = data[3]
                SR = data[4]
                AR = data[5]
                PR = data[6]
            continue
        if data[0] == 'IT':
            wl.append(data[1])
            IT.append(data[2])
            continue
        if data[0] == 'GC':
            data.remove('GC')
            comments = data
if glass_count!=0:WriteYML()


