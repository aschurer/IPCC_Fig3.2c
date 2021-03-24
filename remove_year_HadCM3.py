import iris

HadCM3=iris.load('/exports/csce/datastore/geos/users/aschurer/IPCC_past1000_figure/CMIP5/RawData/tas_Amon_HadCM3_past1000_r1i1p1_185001-200011.nc')[0]
HadCM3=HadCM3[12:]
iris.save(HadCM3, '/exports/csce/datastore/geos/users/aschurer/IPCC_past1000_figure/CMIP5/RawData/tas_Amon_HadCM3_past1000_r1i1p1_185101-200011.nc')
