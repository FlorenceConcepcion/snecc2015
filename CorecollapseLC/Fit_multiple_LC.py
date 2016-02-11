import sncosmo
import astropy
import os
import numpy as np
import matplotlib.pyplot as plt

'''
Note, the following are a small sample of the available built in models in sncosmo:
TypeIb model: source=sncosmo.get_source('snana-2006ep', version='1.0')
TypeIb/c model: sncosmo.get_source('nugent-sn1bc', version='1.1')
TypeIc model: sncosmo.get_source('snana-2004fe',version='1.0')
TypeIIP model: sncosmo.get_source('s11-2004hx',version='1.0')
'''

# Registering ptf48r filter
import os
import sys

sys.path.append(os.path.abspath('..'))
import Bandpass

'''
# Creating a dictionary from a list of Supernovae with their name and type of snhml = open('List/SNlist.txt')
lines = [line.split(' ') for line in hml if len(line) > 2]
types = {}
for line in lines:
    for i in range(len(line)):
        if line[i] == 'SN':
            types[line[0]] = line[i + 1].strip(',')
'''
# Dictionary created from given list of core collapse supernova (using code above)
types = {'10hgi': 'Ib/c', '10pjg': 'IIP', '10qwz': 'IIP', '10rin': 'IIP', '10wmf': 'IIP', '10hdq': 'IIP',
         '10ttd': 'IIP', '10qqd': 'Ic', '10qts': 'Ic-BL', '10tpa': 'IIP', '10xfh': 'Ic', '10osr': 'IIP', '10pbi': 'Ib',
         '10wal': 'Ic', '10vgv': 'Ic-BL', '10hie': 'Ic', '10rmn': 'IIP', '10xik': 'Ic', '10nbf': 'IIP', '10glp': 'IIP',
         '10xem': 'Ic-BL', '10lbo': 'Ic', '10yow': 'Ic', '10osn': 'Ic', '10rem': 'IIP', '10hny': 'IIP', '10hpa': 'IIP',
         '10qif': 'Ib', '10tqv': 'Ic', '10vdl': 'IIP', '10uqn': 'IIP', '10inj': 'Ib', '10rjs': 'IIP', '10npd': 'IIP',
         '10hfe': 'Ic', '10qob': 'IIP', '10ood': 'Ic', '10uhf': 'Ic', '10svt': 'Ic', '10hyq': 'IIP', '10vwg': 'Ic',
         '10tqi': 'Ic', '10uhn': 'Ic', '10myz': 'IIP', '10wve': 'IIP', '10ysd': 'Ic-BL', '10vnv': 'Ib', '10hcw': 'Ib',
         '10tff': 'IIP', '10gtn': 'IIP', '10gvb': 'Ic-BL', '10kui': 'Ib', '10xjv': 'IIP', '10xjr': 'Ib',
         '10iue': 'Ic-BL'}


def fit_curve(hml_dat, source, hml):
    # Retrieving model from source, (requires internet access)
    model = sncosmo.Model(source=source)
    # Fitting model: model parameter z bound.
    res, fitted_model = sncosmo.fit_lc(hml_dat, model, ['z', 't0', 'amplitude'], bounds={'z': (0.005, 0.35)})

    # The following excludes data sets with fewer than 4 data points and not enough distribution..

    # print 'peak', float(res.parameters[1])
    j = 0
    k = 0
    for i in hml[:, 1]:
        if float(i) > float(res.parameters[1]):
            j += 1
        else:
            k += 1

    if j >= 2 and k >= 2:
        print len(hml[:, 1])
        if len(hml[:, 1]) > 3:
            # print res.errors
            # print model

            sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(
                hml[:, 0][0] + ' Type' + types[hml[:, 0][0]] + '\n' + 'Model name: ' + source.name), xfigsize=10)
            plt.show()
            print 'Done:', hml[:, 0][0], 'z:', float(
                res.parameters[0]), 'Reduced chi^2:', res.chisq / res.ndof, 'Data points:', len(
                hml[:, 1]), '\n'  # 'Dof:', res.ndof

        else:
            pass


def select_model(x):
    # Load file containing data set
    hml = np.load(os.path.abspath('') + '/ccdata/' + x)
    # print len(hml[:,1]), 'data points'             #number of rows and therefore data points

    # Adding zpsys and filter columns. zp system used is ab and filter used is ptf48r
    ab = np.zeros(len(hml[:, 1]), dtype='|S2')
    band = np.zeros(len(hml[:, 1]), dtype='|S6')
    for i in range(len(ab)):
        ab[i] = 'ab'
        band[i] = 'ptf48r'
    hml = np.column_stack((hml, ab, band))

    # Converting into a table using astropy with titles: ptfname, time,
    # magnitude, mag_err, flux, flux_err, zp_new, zp, ra, dec, zpsys and filter
    hml_dat = astropy.table.Table(data=hml, names=(
            'ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys',
            'filter'),
                                      dtype=(
                                          'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',
                                            'float', 'str', 'str'))
    #sn types are Ib, Ib/c, Ic, Ic-BL, IIP
    sn_type = types[hml[:, 0][0]]

    # Selecting only Type II supernovae
    if sn_type[1] == 'I':
        # Selecting model
        source = sncosmo.get_source('s11-2004hx',version='1.0')
        fit_curve(hml_dat, source, hml)

    # Selecting only Type Ic supernovae
    elif sn_type[1]== 'c':
        # Selecting model
        source = sncosmo.get_source('snana-2004fe',version='1.0')
        fit_curve(hml_dat, source, hml)

    # Selecting TypeIb sn only
    elif len(sn_type)== 2 and (sn_type[1]!= 'b'):
        # Selecting model
        source = sncosmo.get_source('snana-2006ep', version='1.0')
        fit_curve(hml_dat, source, hml)

    # Selecting TypeIbc sn only
    elif len(sn_type)> 2 and (sn_type[1]!= 'b'):
        # Selecting model
        source = sncosmo.get_source('nugent-sn1bc', version='1.1')
        fit_curve(hml_dat, source, hml)

'''
# Lists all supernova data sets available in CorecollapseLC/ccdata/
cclist = []
for i in os.walk('ccdata/'):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)
print len(cclist), 'supernovae in list'
'''
cclist_retrieved = ['10glp_LC.npy', '10gtn_LC.npy', '10gvb_LC.npy', '10hcw_LC.npy', '10hdq_LC.npy', '10hfe_LC.npy',
                    '10hgi_LC.npy', '10hie_LC.npy', '10hny_LC.npy', '10hpa_LC.npy', '10hyq_LC.npy', '10inj_LC.npy',
                    '10iue_LC.npy', '10kui_LC.npy', '10lbo_LC.npy', '10myz_LC.npy', '10nbf_LC.npy', '10npd_LC.npy',
                    '10ood_LC.npy', '10osn_LC.npy', '10osr_LC.npy', '10pbi_LC.npy', '10pjg_LC.npy', '10qif_LC.npy',
                    '10qob_LC.npy', '10qqd_LC.npy', '10qts_LC.npy', '10qwz_LC.npy', '10rem_LC.npy', '10rin_LC.npy',
                    '10rjs_LC.npy', '10rmn_LC.npy', '10svt_LC.npy', '10tff_LC.npy', '10tpa_LC.npy', '10tqi_LC.npy',
                    '10tqv_LC.npy', '10ttd_LC.npy', '10uhf_LC.npy', '10uhn_LC.npy', '10uqn_LC.npy', '10vdl_LC.npy',
                    '10vgv_LC.npy', '10vnv_LC.npy', '10vwg_LC.npy', '10wal_LC.npy', '10wmf_LC.npy', '10wve_LC.npy',
                    '10xem_LC.npy', '10xfh_LC.npy', '10xik_LC.npy', '10xjr_LC.npy', '10xjv_LC.npy', '10yow_LC.npy',
                    '10ysd_LC.npy']
cclist = cclist_retrieved
for i in cclist[:]:
    select_model(i)
