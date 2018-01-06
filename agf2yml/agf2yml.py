#!/usr/bin/python
# coding: utf-8
###############################################################################
#         THIS PROGRAM IS IN PUBLIC DOMAIN                                    #
#         COPYRIGHT AND RELATED RIGHTS WAIVED VIA CC0 1.0                     #
###############################################################################

#             agf2yml (Zemax to refractiveindex.info converter)

#------------------------------------------------------------------------------
#  dependencies: Python3, NumPy
#------------------------------------------------------------------------------


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CHECK AGF FILE ENCODING !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import os, math

agffiles = []
ymldirs = []
refs = []

agf_catalogs = {}

agf_catalogs['schott'] = {'file': 'input/schott_2017-01-20b.agf',
                          'dir': 'output/schott',
                          'refs': '<a href=\\"http://refractiveindex.info/download/data/2017/schott_2017-01-20b.agf\\">SCHOTT Zemax catalog 2017-01-20b</a> (obtained from <a href=\\"http://www.schott.com/advanced_optics/english/download/\\">http://www.schott.com</a>)<br>See also <a href=\\"http://refractiveindex.info/download/data/2017/schott_2017-01-20.pdf\\">SCHOTT glass data sheets</a>'}

agf_catalogs['ohara'] = {'file': 'input/ohara_2017-11-30.agf',
                         'dir': 'output/ohara', 
                         'refs': '<a href=\\"http://refractiveindex.info/download/data/2017/ohara_2017-11-30.agf\\">OHARA Zemax catalog 2017-11-30</a> (obtained from <a href=\\"http://www.ohara-inc.co.jp/en/product/optical/list/\\">http://www.ohara-inc.co.jp</a>)<br>See also <a href=\\"http://refractiveindex.info/download/data/2017/ohara_2017-11-30.pdf\\">OHARA glass data sheets</a>'}

agf_catalogs['hikari'] = {'file': 'input/nikon_2017-11.agf',
                          'dir': 'output/hikari',
                          'refs': '<a href=\\"http://refractiveindex.info/download/data/2017/nikon_2017-11.agf\\">NIKON Zemax catalog 2017-11</a> (obtained from <a href=\\"http://www.nikon.com/products/glass/lineup/materials/optical/\\">http://www.nikon.com</a>)<br>See also <a href=\\"http://refractiveindex.info/download/data/2017/hikari_2017-04-01.pdf\\">HIKARI glass data sheets</a>'}

agf_catalogs['hoya'] = {'file': 'input/hoya_2017-04-01.agf',
                        'dir': 'output/hoya',
                        'refs': '<a href=\\"http://refractiveindex.info/download/data/2017/hoya_2017-04-01.agf\\">HOYA Zemax catalog 2017-04-01</a> (obtained from <a href=\\"http://www.hoya-opticalworld.com/english/datadownload/\\">http://www.hoya-opticalworld.com</a>)'}

agf_catalogs['sumita'] = {'file': 'input/sumita_2017-02-02.agf', 
                          'dir': 'output/sumita',
                          'refs': '<a href=\\"http://refractiveindex.info/download/data/2017/sumita_2017-02-02.agf\\">SUMITA Zemax catalog 2017-02-02</a> (obtained from <a href=\\"http://www.sumita-opt.co.jp/en/download/\\">http://www.sumita-opt.co.jp</a>)<br>See also <a href=\\"http://refractiveindex.info/download/data/2016/sumita_2016-02-01.pdf\\">SUMITA optical glass data book</a>'}

agf_catalogs['cdgm'] = {'file': 'input/cdgm_2017-09.agf', 
                          'dir': 'output/cdgm',
                          'refs': '<a href=\\"http://refractiveindex.info/download/data/2017/cdgm_2017-09.agf\\">CDGM Zemax catalog 2017-09</a> (obtained from <a href=\\"http://www.cdgmgd.com/208/\\">http://www.cdgmgd.com</a>)'}

class GlassData:
    wl = None
    IT = None
    thickness = None
    
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
    def __init__(self):
        self.wl, self.IT, self.thickness = [],[],[]


def WriteYML(gd, ymldir, references, glass_count): 
    print('{}: {}'.format(glass_count, gd.name))
    yml_file_path = os.path.join(ymldir, gd.name.replace('*','star'))
    try:
        ymlfile = open('{}.yml'.format(yml_file_path), 'w+', encoding='utf-8')
    except TypeError: # Python2 has no encoding= argument
        ymlfile = open('{}.yml'.format(yml_file_path), 'w+')

    ymlfile.write("# this file is part of refractiveindex.info database\n# refractiveindex.info database is in the public domain\n# copyright and related rights waived via CC0 1.0\n\n")
    ymlfile.write('REFERENCES: "' + references + '"\n')
    
    if gd.comments:
        ymlfile.write('COMMENTS: "' + ' '.join(gd.comments) + '"\n')
    
    ymlfile.write('DATA:\n')

    if gd.formula == "1":
        ymlfile.write('  - type: formula 3 \n')
        ymlfile.write('    range: {} {}\n'.format(float(gd.wlmin), float(gd.wlmax)))
        ymlfile.write('    coefficients:')
        coeff_list = [None, 2, -2, -4, -6, -8]
        n_coeffs = len(gd.disp_formula_coefficients)
        for i, k in zip(range(n_coeffs), coeff_list[:n_coeffs]):
            if float(gd.disp_formula_coefficients[i]):
                ymlfile.write(' {}'.format(float(gd.disp_formula_coefficients[i]), k))
                if k != None:
                    ymlfile.write(' {}'.format(k))
                    
    elif gd.formula == "2" or "6":
        ymlfile.write('  - type: formula 2 \n')
        ymlfile.write('    range: {} {}\n'.format(float(gd.wlmin), float(gd.wlmax)))
        ymlfile.write('    coefficients:')
        ymlfile.write(' 0')
        n_coeff = len(gd.disp_formula_coefficients) 
        for i in range(0, 16, 2):
            if n_coeff > i+1 and float(gd.disp_formula_coefficients[i]):
                ymlfile.write(' {} {}'.format(float(gd.disp_formula_coefficients[i]),
                                              float(gd.disp_formula_coefficients[i+1])))
                
    elif gd.formula == "3":
        ymlfile.write('  - type: formula 4 \n')
        ymlfile.write('    range: {} {}\n'.format(float(gd.wlmin), float(gd.wlmax)))
        ymlfile.write('    coefficients:')
        n_coeffs = len(gd.disp_formula_coefficients)
        ymlfile.write(' {} {} 2 {} 2 {} 0 {} 2'.format(float(gd.disp_formula_coefficients[0])+1,
                                                       float(gd.disp_formula_coefficients[1]),
                                                       float(gd.disp_formula_coefficients[2]),
                                                       float(gd.disp_formula_coefficients[3]),
                                                       float(gd.disp_formula_coefficients[4])))
    elif gd.formula == "11":
        ymlfile.write('  - type: formula 3 \n')
        ymlfile.write('    range: {} {}\n'.format(float(gd.wlmin), float(gd.wlmax)))
        ymlfile.write('    coefficients:')
        coeff_list = [None, 2, -2, -4, -6, -8, -10, -12]
        n_coeffs = len(gd.disp_formula_coefficients)
        for i, k in zip(range(n_coeffs), coeff_list[:n_coeffs]):
            if float(gd.disp_formula_coefficients[i]):
                ymlfile.write(' {}'.format(float(gd.disp_formula_coefficients[i]), k))
                if k != None:
                    ymlfile.write(' {}'.format(k))
    
    elif gd.formula == "12":
        ymlfile.write('  - type: formula 3 \n')
        ymlfile.write('    range: {} {}\n'.format(float(gd.wlmin), float(gd.wlmax)))
        ymlfile.write('    coefficients:')
        coeff_list = [None, 2, -2, -4, -6, -8, 4, 6]
        n_coeffs = len(gd.disp_formula_coefficients)
        for i, k in zip(range(n_coeffs), coeff_list[:n_coeffs]):
            if float(gd.disp_formula_coefficients[i]):
                ymlfile.write(' {}'.format(float(gd.disp_formula_coefficients[i]), k))
                if k != None:
                    ymlfile.write(' {}'.format(k))
    
    elif gd.formula == "13":
        ymlfile.write('  - type: formula 3 \n')
        ymlfile.write('    range: {} {}\n'.format(float(gd.wlmin), float(gd.wlmax)))
        ymlfile.write('    coefficients:')
        coeff_list = [None, 2, 4, -2, -4, -6, -8, -10, -12]
        n_coeffs = len(gd.disp_formula_coefficients)
        for i, k in zip(range(n_coeffs), coeff_list[:n_coeffs]):
            if float(gd.disp_formula_coefficients[i]):
                ymlfile.write(' {}'.format(float(gd.disp_formula_coefficients[i]), k))
                if k != None:
                    ymlfile.write(' {}'.format(k))
                
    ymlfile.write('\n')
    
    ymlfile.write('  - type: tabulated k\n')
    ymlfile.write('    data: |\n')  
    for i in range (0, len(gd.wl)):
        if float(gd.IT[i]) != 0:
            k = -float(gd.wl[i]) / (4*math.pi) * \
                math.log(float(gd.IT[i])) / (float(gd.thickness[i])*1000)
            ymlfile.write('        {:.3f} {:.4E}\n'.format(float(gd.wl[i]), k))
    
    ymlfile.write('SPECS:\n')
    ymlfile.write('    n_is_absolute: false\n')
    ymlfile.write('    λ_is_vacuum: false\n')
    if len(gd.thermal_disp_coefficients) > 6 and gd.thermal_disp_coefficients[6] != '':
        ymlfile.write('    temperature: {:.1f} °C\n'.format(float(gd.thermal_disp_coefficients[6])))
    if len(gd.thermal_disp_coefficients) > 5 and \
            gd.thermal_disp_coefficients[0] != '' and \
            gd.thermal_disp_coefficients[1] != '' and \
            gd.thermal_disp_coefficients[2] != '' and \
            gd.thermal_disp_coefficients[3] != '' and \
            gd.thermal_disp_coefficients[4] != '' and \
            gd.thermal_disp_coefficients[5] != '' and \
            (float(gd.thermal_disp_coefficients[0]) or
             float(gd.thermal_disp_coefficients[1]) or
             float(gd.thermal_disp_coefficients[2]) or
             float(gd.thermal_disp_coefficients[3]) or
             float(gd.thermal_disp_coefficients[4]) or
             float(gd.thermal_disp_coefficients[5])): 
        ymlfile.write('    coefficients_of_thermal_dispersion: {}\n'\
            .format(' '.join([format(float(gd.thermal_disp_coefficients[i])) \
                              for i in range(6) ]))) 
    ymlfile.write('    n_d: {}\n'.format(gd.nd))
    ymlfile.write('    V_d: {}\n'.format(gd.vd))
    if float(gd.glasscode) > 100000:
        ymlfile.write('    glass_code: {}\n'.format(gd.glasscode))
    
    status_name = {1:'standard', 2:'preferred', 3:'special', 4:'obsolete',
                   5:'melt'}
    for num, desc in status_name.items():
        if int(float(gd.status)) == num:
             ymlfile.write('    glass_status: {}\n'.format(desc))

    if gd.meltfreq != '' and gd.meltfreq != 0:
        ymlfile.write('    glass_melt_frequency: {}\n'.format(gd.meltfreq))
    if gd.density != '' and gd.density != '-' and float(gd.density):
        ymlfile.write('    density: {} g/cm<sup>3</sup>\n'.format(float(gd.density)))
    if gd.cte1 != '' and gd.cte1 != '-' and float(gd.cte1): #-30...70 C
        ymlfile.write('    coefficient_of_thermal_expansion_1: {} K<sup>-1</sup>\n'\
                      .format(round(float(gd.cte1)*1.0e-6,15)))
    if gd.cte2 != '' and gd.cte2 != '-' and float(gd.cte2): #20...300 C
        ymlfile.write('    coefficient_of_thermal_expansion_2: {} K<sup>-1</sup>\n'.\
                       format(round(float(gd.cte2)*1.0e-6,15)))
    if gd.dpgf != '' and gd.dpgf != '-' and float(gd.dpgf):
        ymlfile.write('    ΔP_gF: {}\n'.format(float(gd.dpgf)))
    if gd.CR != '' and gd.CR != '-' and float(gd.CR) != -1:
        ymlfile.write('    climatic_resistance: {}\n'.format(float(gd.CR)))
    if gd.FR != '' and gd.FR != '-' and float(gd.FR) != -1:
        ymlfile.write('    stain_resistance: {}\n'.format(float(gd.FR)))
    if gd.SR != '' and gd.SR != '-' and float(gd.SR) != -1:
        ymlfile.write('    acid_resistance: {}\n'.format(float(gd.SR)))
    if gd.AR != '' and gd.AR != '-' and float(gd.AR) != -1:
        ymlfile.write('    alkali_resistance: {}\n'.format(float(gd.AR)))
    if gd.PR != '' and gd.PR != '-' and float(gd.PR) != -1:
        ymlfile.write('    phosphate_resistance: {}\n'.format(float(gd.PR)))

    ymlfile.close()
    
    print('ok')

def process(agf_file, out_dir, ref):
    glass_count = 0
    gd = GlassData()
    with open(agf_file) as agf:
        for line in agf:
            data = [k.strip() for k in line.split(" ")]
            if not len(data):
                continue
            if data[0] == 'NM':
                if glass_count!=0:
                    WriteYML(gd, out_dir, ref, glass_count)
                    #sys.exit(0)
                glass_count +=1
                gd = GlassData()
                gd.name = data[1]
                
                gd.formula = data[2]
                gd.glasscode = data[3]
                gd.nd = data[4]
                gd.vd = data[5]
                if len(data) >= 8: gd.status = data[7] 
                if len(data) >= 9 and data.count('-') < 0: gd.meltfreq = data[8]
                else: gd.meltfreq = 0
                continue
            if data[0] == 'LD':
                gd.wlmin = data[1]
                gd.wlmax = data[2]
                continue
            if data[0] == 'CD':
                gd.disp_formula_coefficients = data[1:]
                continue
            if data[0] == 'TD':
                gd.thermal_disp_coefficients = data[1:]
                continue
            if data[0] == 'ED':
                gd.cte1 = data[1]
                gd.cte2 = data[2]
                gd.density = data[3]
                gd.dpgf = data[4]
                continue
            if data[0] == 'OD':
                #? = data[1]
                if len(data)>2:
                    gd.CR = data[2]
                    if gd.CR == '3-4':
                        gd.CR = '3.5'
                    gd.FR = data[3]
                    gd.SR = data[4]
                    gd.AR = data[5]
                    gd.PR = data[6]
                continue
            if data[0] == 'IT' and data[1]!='' and data[2]!='' and data[3]!='':
                gd.wl.append(data[1])
                gd.IT.append(data[2])
                gd.thickness.append(data[3])
                continue
            if data[0] == 'GC':
                gd.comments = data[1:]

        if glass_count!=0:
            WriteYML(gd, out_dir, ref, glass_count)

### main program
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python {} <catalog name>".format(__file__))
        print("Catalog name can be on of the following: ")
        for name in agf_catalogs.keys():
            print("{} ".format(name))
        sys.exit(1)
    catalog = agf_catalogs[sys.argv[1]]
    # Make output dir
    if not os.path.exists(catalog['dir']):
        os.mkdir(catalog['dir'])
    process(catalog['file'], catalog['dir'], catalog['refs'])
    