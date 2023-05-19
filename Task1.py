import pyspedas
import pytplot

# ---- TASK 1 --------
# TASK 1.1

trange=['2018-11-11/06:00', '2018-11-12/12:00']

fields_vars = pyspedas.psp.fields(trange=trange, datatype='mag_rtn_1min', level='l2', time_clip=True)

# ------ MAGNETIC FIELD --------

Bsep=pytplot.split_vec('psp_fld_l2_mag_RTN_1min')
B_r=pytplot.get_data('psp_fld_l2_mag_RTN_1min_x')
B_t=pytplot.get_data('psp_fld_l2_mag_RTN_1min_y')
B_n=pytplot.get_data('psp_fld_l2_mag_RTN_1min_z')
B=(B_r.y**2+B_t.y**2+B_n.y**2)**(1/2) 
pytplot.store_data('Bt', data={'x': B_r.times, 'y': B})

pytplot.join_vec(['psp_fld_l2_mag_RTN_1min_x','psp_fld_l2_mag_RTN_1min_y','psp_fld_l2_mag_RTN_1min_z','Bt'],new_tvar='b_mult')

# ---- SOLAR WIND VELOCITY (MOMENT FITS) ------

# Do they just use SPC l3 proton moments?

# Moments (it's this one)

spc_vars = pyspedas.psp.spc(trange=trange, datatype='l3i', level='l3', time_clip=True)
pytplot.split_vec('psp_spc_vp_moment_RTN')

# Fits averaged (alphas seem to give NaN)

n_p=pytplot.get_data('psp_spc_np_fit')
n_a=pytplot.get_data('psp_spc_na_fit')
n_3=pytplot.get_data('psp_spc_n3_fit')

v_p=pytplot.get_data('psp_spc_vp_fit_RTN')
pytplot.split_vec('psp_spc_vp_fit_RTN')
v_pr=pytplot.get_data('psp_spc_vp_fit_RTN_x')

v_a=pytplot.get_data('psp_spc_va_fit_RTN')
pytplot.split_vec('psp_spc_va_fit_RTN')
v_ar=pytplot.get_data('psp_spc_va_fit_RTN_x')

v_3=pytplot.get_data('psp_spc_v3_fit_RTN')
pytplot.split_vec('psp_spc_v3_fit_RTN')
v_3r=pytplot.get_data('psp_spc_v3_fit_RTN_x')

pytplot.store_data('Vw_avgd', data={'x': v_p.times, 'y': (v_3r.y*n_3.y+v_ar.y*n_a.y+v_pr.y*n_p.y)/(n_3.y+n_a.y+n_p.y)})

# ------ SOLAR WIND DENSITY -----

# Already obtained for the fits attempt

# ----- PROTON TEMPERATURE ------

mp=1.6726E-27
k=1.380649E-23
spi_vars_p = pyspedas.psp.spi(trange=trange, datatype='spi_sf00_l3_mom', level='l3', time_clip=True)
wp=pytplot.get_data('psp_spc_wp_moment')
t=(wp.y**2*mp)/(2*k)
pytplot.store_data('T_p', data={'x': wp.times, 'y': t})

# ----- BETA VALUE ------

Bsq=B**2
pytplot.store_data('Bsq',data={'x': B_r.times, 'y': Bsq})
pytplot.tplot_math.multiply('psp_spc_np_fit','T_p','npxt')
pytplot.tplot_math.divide('npxt','Bsq','npxt/Bsq')
npxtBsq=pytplot.get_data('npxt/Bsq')


pytplot.store_data('beta',data={'x': npxtBsq.times, 'y': npxtBsq.y*3.5*10})
pytplot.options('beta', 'yrange', [0,3])

# ----- FINAL PLOTS -----

# TASK 1.2

pytplot.timebar(1541980260, color='green')
pytplot.timebar(1542003420, color='green')
pytplot.timebar(1541959200, color='red')

pytplot.options('beta', 'yrange', [0,3])
pytplot.options('psp_spc_np_fit', 'yrange', [0,600])
pytplot.tplot(['b_mult','Bt','psp_spc_vp_moment_RTN_x','psp_spc_np_fit','T_p','beta'])

# TASK 1.3

# Propagation speed (+density)

trange=['2018-11-11/23:51','2018-11-12/06:17']
spc_vars = pyspedas.psp.spc(trange=trange, datatype='l3i', level='l3', time_clip=True)
pytplot.split_vec('psp_spc_vp_moment_RTN')
vp=pytplot.get_data('psp_spc_vp_moment_RTN')
pytplot.tplot_math.avg_res_data('psp_spc_vp_moment_RTN_x', len(vp.y), 'vprop')
vprop=pytplot.get_data('vprop')

# Density ratio

# MO density

np_mo=pytplot.get_data('psp_spc_np_fit')
pytplot.tplot_math.avg_res_data('psp_spc_np_fit', len(np_mo.y), 'nmo')
nmo=pytplot.get_data('nmo')


# Sw density

trange=['2018-11-10/18:00','2018-11-11/18:00']
spc_vars = pyspedas.psp.spc(trange=trange, datatype='l3i', level='l3', time_clip=True)
np_sw=pytplot.get_data('psp_spc_np_fit')
pytplot.tplot_math.avg_res_data('psp_spc_np_fit', len(np_sw.y), 'nsw')
nsw=pytplot.get_data('nsw')

# TASK 1.4

# How radial distance varies during the MO

trange=['2018-11-11/23:51','2018-11-12/06:17']
spc_vars = pyspedas.psp.spc(trange=trange, datatype='l3i', level='l3', varnames=['sc_pos_HCI','sc_vel_HCI','vp_moment_RTN'], time_clip=True)
pytplot.tplot_rename('psp_spc_sc_pos_HCI','dsc')
pytplot.split_vec('dsc')

d_x = pytplot.get_data('dsc_x')
d_y = pytplot.get_data('dsc_y')
d_z = pytplot.get_data('dsc_z')
d = (d_x.y**2+d_y.y**2+d_z.y**2)**(1/2)

pytplot.store_data('d', data={'x': d_x.times, 'y': d})

pytplot.tplot_math.avg_res_data('d', len(d), 'dav')
dav=pytplot.get_data('dav')

# TASK 1.5

# Comparing flow velocity with PSP velocity

pytplot.tplot_rename('psp_spc_sc_vel_HCI','vsc')
pytplot.split_vec('vsc')

v_xsc = pytplot.get_data('vsc_x')
v_ysc = pytplot.get_data('vsc_y')
v_zsc = pytplot.get_data('vsc_z')
vsc = (v_xsc.y**2+v_ysc.y**2+v_zsc.y**2)**(1/2)

pytplot.store_data('vsctot', data={'x': v_xsc.times, 'y': vsc})

pytplot.split_vec('psp_spc_vp_moment_RTN')

v_x = pytplot.get_data('psp_spc_vp_moment_RTN_x')
v_y = pytplot.get_data('psp_spc_vp_moment_RTN_y')
v_z = pytplot.get_data('psp_spc_vp_moment_RTN_z')
v = (v_x.y**2+v_y.y**2+v_z.y**2)**(1/2)

pytplot.store_data('vtot', data={'x': v_x.times, 'y': v})

pytplot.join_vec(['vsctot','vtot'],new_tvar='vtotcomparison')

pytplot.tplot_math.divide('vtot','vsctot','vratio')

pytplot.tplot(['vratio'])

