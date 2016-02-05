#!/usr/bin/python
# coding: utf-8
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

agffile = 'input/schottzemax-20150722.agf'
ymldir = 'output/schott'
references = '1) <a href=\\"http://refractiveindex.info/download/data/2015/schott-optical-glass-collection-datasheets-july-2015-us.pdf\\">SCHOTT optical glass data sheets 2015-07-22</a><br>2) <a href=\\"http://refractiveindex.info/download/data/2015/schottzemax-20150722.agf\\">SCHOTT Zemax catalog 2015-07-22</a>'

#agffile = 'input/OHARA_151201.agf'
#ymldir = 'output/ohara'
#references = '1) <a href=\\"http://refractiveindex.info/download/data/2015/ohara_2015-12-01.pdf\\">OHARA optical glass datasheets 2015-12-01</a><br>2) <a href=\\"http://refractiveindex.info/download/data/2015/OHARA_151201.agf\\">OHARA Zemax catalog 2015-12-01</a>'

#agffile = 'input/HIKARI.agf'
#ymldir = 'output/hikari'
#references = '1) <a href=\\"http://refractiveindex.info/download/data/2015/HIKARI_Catalog.pdf\\">HIKARI optical glass catalog 2015-04-01</a><br>2) <a href=\\"http://refractiveindex.info/download/data/2015/HIKARI.agf\\">HIKARI Zemax catalog</a>'

#agffile = 'input/HOYA20150618.agf'
#ymldir = 'output/hoya'
#references = '<a href=\\"http://refractiveindex.info/download/data/2015/HOYA20150618.agf\\">HOYA Zemax catalog 2015-06-18</a>'

wl = []
IT = []
thickness = []


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
    del thickness[:]



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
    ymlfile.write('    range: {} {}\n'.format(wlmin, wlmax))
    ymlfile.write('    coefficients:')
    if formula == "1" or formula == "13":
        for i, k in zip(range(8), ['', 2, 4, -2, -6, -8, -10, -12]):
            if float(disp_formula_coefficients[i]):
                ymlfile.write(' {} {}'.format(disp_formula_coefficients[i], k))
    if formula == "2":
        ymlfile.write(' 0')
        num_coeff = len(disp_formula_coefficients) 
        for i in range(0, 16, 2):
            if num_coeff > i+1 and float(disp_formula_coefficients[i]):
                ymlfile.write(' {} {}'.format(disp_formula_coefficients[i],
                                              disp_formula_coefficients[i+1]))
    ymlfile.write('\n')
    
    ymlfile.write('  - type: tabulated k\n')
    ymlfile.write('    data: |\n')  
    for i in range (0,len(wl)):
        if float(IT[i]) != 0:
            k = -float(wl[i]) / (4*math.pi) * \
                math.log(float(IT[i])) / (float(thickness[i])*1000)
            ymlfile.write('        {:.3f} {:.4E}\n'.format(float(wl[i]), k))
    
    ymlfile.write('INFO:\n')
    ymlfile.write('    n_is_absolute: false\n')
    ymlfile.write('    λ_is_vacuum: false\n')
    if len(thermal_disp_coefficients) > 6:
        ymlfile.write('    temperature: {}  °C\n'.format(thermal_disp_coefficients[6]))
    if len(thermal_disp_coefficients) > 5 and \
        (float(thermal_disp_coefficients[0]) or \
            float(thermal_disp_coefficients[1]) or \
            float(thermal_disp_coefficients[2]) or \
            float(thermal_disp_coefficients[3]) or \
            float(thermal_disp_coefficients[4]) or \
            float(thermal_disp_coefficients[5])): 
        ymlfile.write('    coefficients_of_thermal_dispersion: {}\n'\
            .format(' '.join([format(float(thermal_disp_coefficients[i])) \
                              for i in range(6) ]))) 
    ymlfile.write('    n_d: {}\n'.format(nd))
    ymlfile.write('    V_d: {}\n'.format(vd))
    if float(glasscode) > 100000:
        ymlfile.write('    glass_code: ' + glasscode + '\n')
    
    status_name = {1:'standard', 2:'preferred', 3:'special', 4:'obsolete',
                   5:'melt'}
    for num, desc in status_name.items():
        if int(float(status)) == num:
             ymlfile.write('    glass_status: {}\n'.format(desc))

    if meltfreq != 0:
        ymlfile.write('    glass_melt_frequency: {}\n'.format(meltfreq))
    if density != '-' and float(density):
        ymlfile.write('    density: {} g/cm<sup>3</sup>\n'.format(density))
    if cte1 != '-' and float(cte1): #-30...70 C
        ymlfile.write('    coefficient_of_thermal_expansion_1: {} K<sup>-1</sup>\n'\
                      .format(round(float(cte1)*1.0e-6,15)))
    if cte2 != '-' and float(cte2): #20...300 C
        ymlfile.write('    coefficient_of_thermal_expansion_2: {} K<sup>-1</sup>\n'.\
                       format(round(float(cte2)*1.0e-6,15)))
    if dpgf != '-' and float(dpgf):
        ymlfile.write('    ΔP_gF: {}\n'.format(dpgf))
    if CR != '-' and float(CR) != -1:
        ymlfile.write('    climatic_resistance: {}\n'.format(CR))
    if FR != '-' and float(FR) != -1:
        ymlfile.write('    stain_resistance: {}\n'.format(FR))
    if SR != '-' and float(SR) != -1:
        ymlfile.write('    acid_resistance: {}\n'.format(float(SR)))
    if AR != '-' and float(AR) != -1:
        ymlfile.write('    alkali_resistance: {}\n'.format(AR))
    if PR != '-' and float(PR) != -1:
        ymlfile.write('    phosphate_resistance: {}\n'.format(PR))

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
            thickness.append(data[3])
            continue
        if data[0] == 'GC':
            data.remove('GC')
            comments = data
if glass_count!=0:WriteYML()


