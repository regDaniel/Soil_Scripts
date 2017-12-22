load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; order of loading is important
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------------
begin
   case="7_years"

   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+case+"/output/out01/1979/"

   in     = addfile(folder+"out01s.nc","r")
   ;s      = in->S_SO ; saturation
   wso      = in->W_SO ; saturation
   wio      = in->W_SO_ICE ; soil ice
   wii      = in->W_I ; intercept
   wsnow    = in->W_SNOW ; snow
   ;klw    = in->KLW ; conductivity at layer k-0.5
   ;dlw    = in->DLW ; conductivity at layer k-0.5
;   zh     = in->soil1_bnds(1,:)
;   nz     = dimsizes(zh)
   nz = 10
   level  = ispan(0,nz-1,1)
   ;in     = addfile(folder+"soil_michi.nc","r")
   ;sl     = in->SATLEV ; lowest level that is completely saturated
   ;wt     = in->WTDEPTH ; lowest level that is completely saturated
   ;q      = in->Q  ; flux
   S=dimsizes(wso)
   nx=S(3)
   ny=S(2)

   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+case+"/output/out04/1979/"
   in     = addfile(folder+"out04s.nc","r")
   evap_a  = in->ALHFL_S
   evap_bs = in->ALHFL_BS
   evap_pl = in->ALHFL_PL
   t2m     = in->T_2M

   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+case+"/output/out05/1979/"
   in     = addfile(folder+"out05s.nc","r")
   rs     = -1./6.*in->RUNOFF_S
   rg     = -1./6.*in->RUNOFF_G
   melt   = -1./6.*in->SNOW_MELT
   ;in     = addfile(folder+"soil.nc","r")
   ;melt   = in->SNOW_MELT

   time   = in->time
   timeo  = time/86400.; convert to days
   timeo  = timeo;-365.-366.-15.5; convert to days
   nt     = dimsizes(time)
   print(timeo)
   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+case+"/output/out03/1979/"
   in     = addfile(folder+"out03s_summed.nc","r")
   rain  = (in->RAIN_CON+in->RAIN_GSP)/6. 
   prec  = (in->TOT_PREC) 
   ;raini  = in->RAIN_CON+in->RAIN_GSP 
   if (isfilevar(in, "DEW")) then
      ldew=True
      dew   = in->DEW 
   else
      ldew=False
   end if
   ;preci  = in->TOT_PREC 
   timer  = in->time
   ntr    = dimsizes(timer)

   ;create mask, where we set soil types ice, rock and sea/ sea ice to missing
   folder ="/project/pr04/regenass/cosmo5_validation/"
   in     = addfile(folder+"soiltyp.nc","r")
   soilt  = in->SOILTYP
   maskst = where(((soilt.eq.1).or.(soilt.eq.2).or.(soilt.eq.9).or.(soilt.eq.0)),-999.,1.)
   masks  = maskst(0,:,:)
   masks@_FillValue=-999.
   lat2d     = in->lat
   lon2d     = in->lon
   ;rainsize = dimsizes(raini)
   ;rain = new((/nt,rainsize(1),rainsize(2)/),"float",-999.)
   ;dew  = new((/nt,rainsize(1),rainsize(2)/),"float",-999.)
   ;prec = new((/nt,rainsize(1),rainsize(2)/),"float",-999.)
   ; consider only land points
   ;rain(0,:,:) = raini(0,:,:)
   ;dew(0,:,:)  = dewi(0,:,:)
   ;prec(0,:,:) = preci(0,:,:)
   ;do tt=1,nt-1
   ;   ttr = tt*6
   ;   rain(tt,:,:)   = raini(ttr-5,:,:)+raini(ttr-4,:,:)+raini(ttr-3,:,:)+raini(ttr-2,:,:)+raini(ttr-1,:,:)+raini(ttr,:,:)
   ;;   dew(tt,:,:)   = dewi(ttr-5,:,:)+dewi(ttr-4,:,:)+dewi(ttr-3,:,:)+dewi(ttr-2,:,:)+dewi(ttr-1,:,:)+dewi(ttr,:,:)
   ;;   prec(tt,:,:)   = preci(ttr-5,:,:)+preci(ttr-4,:,:)+preci(ttr-3,:,:)+preci(ttr-2,:,:)+preci(ttr-1,:,:)+preci(ttr,:,:)
   ;end do
   ;rain = rain * 6. ; convert to mm/6h
   
   ; apply mask
   do tt=0,nt-1
      rain(tt,:,:)   = where(ismissing(masks),-999.,rain(tt,:,:))
      if (ldew) then
         dew(tt,:,:)   = where(ismissing(masks),-999.,dew(tt,:,:))
         dew@_FillValue = -999.
      end if
      prec(tt,:,:)   = where(ismissing(masks),-999.,prec(tt,:,:))
   end do
   rain@_FillValue = -999.
   prec@_FillValue = -999.

   do tt=0,nt-1
      evap_a(tt,:,:)   = where(ismissing(masks),-999.,evap_a(tt,:,:))
      evap_bs(tt,:,:)   = where(ismissing(masks),-999.,evap_bs(tt,:,:))
      do kk=0,nz-2
         evap_pl(tt,kk,:,:)   = where(ismissing(masks),-999.,evap_pl(tt,kk,:,:))
         wso(tt,kk,:,:)       = where(ismissing(masks),-999.,wso(tt,kk,:,:))
         wio(tt,kk,:,:)       = where(ismissing(masks),-999.,wio(tt,kk,:,:))
      end do
      wii(tt,:,:)       = where(ismissing(masks),-999.,wii(tt,:,:))
      wsnow(tt,:,:)     = where(ismissing(masks),-999.,wsnow(tt,:,:))
      rg(tt,:,:)        = where(ismissing(masks),-999.,rg(tt,:,:))
      rs(tt,:,:)        = where(ismissing(masks),-999.,rs(tt,:,:))
      melt(tt,:,:)      = where(ismissing(masks),-999.,melt(tt,:,:))
      t2m(tt,:,:)       = where(ismissing(masks),-999.,t2m(tt,:,:))
   end do
   evap_a@_FillValue  = -999.
   evap_bs@_FillValue = -999.
   evap_pl@_FillValue = -999.
   wso@_FillValue     = -999.
   wio@_FillValue     = -999.
   wii@_FillValue     = -999.
   wsnow@_FillValue   = -999.
   rg@_FillValue      = -999.
   rs@_FillValue      = -999.
   melt@_FillValue    = -999.
   t2m@_FillValue     = -999.

   evap_a   = evap_a/2.51e06*3600.;21600. ; change here
   evap_bs  = evap_bs/2.51e06*3600.;21600.; change here
   evap_pl  = evap_pl/2.51e06*3600.;21600.; change here
   rs   = -rs

   dwdt   = new((/nt,nz,ny,nx/),"float",-999)
   dwidt   = new((/nt,ny,nx/),"float",-999)
   dwsdt   = new((/nt,ny,nx/),"float",-999)

   do tt=1,nt-1
      ;dwdt(tt) = dim_sum_n((w(tt,:)-w(tt-1,:))/(time(tt)-time(tt-1)),0)
      dwdt(tt,:,:,:) = (wso(tt,:,:,:)-wso(tt-1,:,:,:))/doubletofloat(time(tt)-time(tt-1))*1000.*3600. ; convert to mm/6h -sat
      dwidt(tt,:,:) = (wii(tt,:,:)-wii(tt-1,:,:))/doubletofloat(time(tt)-time(tt-1))*1000.*3600. ; convert to mm/6h -intercept
      dwsdt(tt,:,:) = (wsnow(tt,:,:)-wsnow(tt-1,:,:))/doubletofloat(time(tt)-time(tt-1))*1000.*3600. ; convert to mm/6h -snow
   end do

   ;dwdt_sum  = dim_sum(dwdt(:,0:7))
   ;w_sum     = dim_sum(wm(:,0:7))
   ;wi_sum    = dim_sum(wim(:,0:7))
   ;empl_sum  = dim_sum(empl(:,0:7))
   ;DR divide everything by 6 to get values in mm/h instead of mm/6h
   dwdt_sum  = -1.*dim_sum_n(dwdt(:,0:8,:,:),1)/6.; change here
   w_sum     = dim_sum_n(wso(:,0:8,:,:),1)/6.; change here
   wi_sum    = dim_sum_n(wii(:,:,:),1)/6.; change here
   empl_sum  = -1.*dim_sum_n(evap_pl(:,0:8,:,:),1); change here

   ;if (ldew) then
   ;   sum_all = pm+rsm+rgm+ema+dm
   ;else
   ;   sum_all = pm+rsm+rgm+ema
   ;end if
   if (ldew) then
      sum_all = (rain+rs+rg+evap_bs+empl_sum+melt+dew)
   else
      sum_all = (rain+rs+rg+evap_bs+empl_sum+melt)
   end if
   ;rel_error = abs(dwdt_sum+dwidt+dwsdt-sum_all)/sum_all*100.
   ;error = dwdt_sum+dwidt+dwsdt-sum_all

   ;rel_error = abs(dwdt_sum-sum_all)/sum_all*100.
   error = dwdt_sum-sum_all

   accum_error = error
   accum_error = 0.
   do tt=1,nt-1
      accum_error(tt,:,:) = accum_error(tt-1,:,:)+error(tt,:,:)
   end do
   ;accum_error = accum_error*4.*30.

   folder ="/project/pr04/regenass/cosmo5_validation/"+case+"/figures/"
   figname = folder+"soil_model_area_"+case
   wks     = gsn_open_wks("ps",figname)

   ; Set resources for plot

   res                     = True

   res@pmLegendDisplayMode = "Always"              ; turn on legend

   res                      = True          ; plot mods desired
   res@gsnFrame             = True          ; False for overlay
   res@gsnDraw              = True          ; False for overlay
   res@gsnMaximize          = True          ; Maximize plot in frame.
   res@gsnSpreadColors      = True

   res@cnFillOn             = True          ; Turn on contour fill
   ;res@cnFillMode           = "RasterFill"  ; Turn on raster fill
   res@cnFillMode           = "CellFill"  ; Turn on raster fill
   res@cnLinesOn            = False         ; Turn off contour lines
   res@cnLineLabelsOn       = False    ; turns off contour line labels
   res@cnLevelSelectionMode = "ManualLevels"
 ;  res@cnMinLevelValF       =  0
 ;  res@cnMaxLevelValF       =  1
 ;  res@cnLevelSpacingF      =  0.1

 ;  res@lbLabelBarOn         = True
   res@lbOrientation        = "vertical"     ; vertical label bar


   res@mpDataBaseVersion = "MediumRes"
   res@mpOutlineBoundarySets = "Geophysical"
   ;res@mpDataBaseVersion     = "Ncarg4_1"
   res@mpDataSetName         = "Earth..4"

   ;res@mpProjection          = "CylindricalEquidistant"
   res@mpProjection          = "Mercator"
   res@mpFillOn              = False
   res@mpGridAndLimbOn       = True             ; turn on lat/lon lines
   res@mpPerimOn             = False            ; turn off box around plot
   res@mpGeophysicalLineThicknessF = 2
   res@mpNationalLineThicknessF    = 2

   res@mpCenterLonF          = 10
   res@mpCenterLatF          = 51
   res@mpLimitMode          = "LatLon"

   res@mpMinLatF            = 24.
   res@mpMaxLatF            = 75.
   res@mpMinLonF            = -12.
   res@mpMaxLonF            = 36.


   res@tiMainString="MOD, bbc2, pompa"

   res@mpMinLatF             = 28.
   res@mpMaxLatF             = 71.
   res@mpMinLonF             = -10.
   res@mpMaxLonF             = 30.

   gsn_define_colormap(wks,"BlueDarkRed18")  ; choose a colormap

   res@cnMinLevelValF       =  -0.4
   res@cnMaxLevelValF       =  0.4
   res@cnLevelSpacingF      =  0.04
   res@gsnLeftString        = "accum error"
   res@gsnRightString       = "m"
   ;delete(res@cnFillColors)

   accum_error@lon2d = lon2d
   accum_error@lat2d = lat2d

   plot = gsn_csm_contour_map(wks,accum_error(nt-1,:,:)/1000.,res)


   xmin=99
   xmax=119

   ymin=57
   ymax=64

   lat2ds=lat2d(ymin:ymax,xmin:xmax)
   lon2ds=lon2d(ymin:ymax,xmin:xmax)

   delete(accum_error@lon2d)
   delete(accum_error@lat2d)
   accum_errors=accum_error(:,ymin:ymax,xmin:xmax)

   accum_errors@lon2d = lon2ds
   accum_errors@lat2d = lat2ds

   plot = gsn_csm_contour_map(wks,accum_errors(nt-1,:,:)/1000.,res)

   rains    = rain(:,ymin:ymax,xmin:xmax)
   rss      = rs(:,ymin:ymax,xmin:xmax)
   rgs      = rg(:,ymin:ymax,xmin:xmax)
   melts    = melt(:,ymin:ymax,xmin:xmax)
   t2ms     = t2m(:,ymin:ymax,xmin:xmax)
   evap_as  = evap_a(:,ymin:ymax,xmin:xmax)
   evap_bss = evap_bs(:,ymin:ymax,xmin:xmax)
   evap_pls = evap_pl(:,:,ymin:ymax,xmin:xmax)

   sum_alls  = sum_all(:,ymin:ymax,xmin:xmax)
   dwdt_sums = dwdt_sum(:,ymin:ymax,xmin:xmax)
   dwidts    = dwidt(:,ymin:ymax,xmin:xmax)

   rm    = dim_avg_n(rains,(/1,2/))
   rsm   = -dim_avg_n(rss,(/1,2/))
   rgm   = -dim_avg_n(rgs,(/1,2/))
   mm    = -dim_avg_n(melts,(/1,2/))
   t2mm  =  dim_avg_n(t2m,(/1,2/))-273.15

   ema   = dim_avg_n(evap_as,(/1,2/))
   embs  = dim_avg_n(evap_bss,(/1,2/))
   empl  = dim_avg_n(evap_pls,(/2,3/))
   empl_sums  = dim_sum(empl(:,0:8)); change here

   sm    = -dim_avg_n(sum_alls,(/1,2/))
   dw    = -dim_avg_n(dwdt_sums,(/1,2/))
   dwi   = -dim_avg_n(dwidts,(/1,2/))
   
   xyres                        = True
   xyres@gsnFrame               = False ; don't advance frame
   xyres@gsnDraw                = False ; don't draw

   xyres@tiXAxisString          = "Time [days]"
   xyres@trYReverse             = True
   xyres@tmXBLabelsOn           = True
   xyres@tmXTLabelsOn           = False
   xyres@tmYRLabelsOn           = False
   xyres@tmYLLabelsOn           = True
   xyres@xyLineThicknessF       = 2
   xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
   xyres@lgLabelPosition        = "Left"
   xyres@pmLegendParallelPosF   = .8                  ; move units right
   xyres@pmLegendWidthF         = 0.08                ; Change width and
   xyres@pmLegendHeightF        = 0.15                ; height of legend.
   xyres@lgPerimOn              = False               ; turn off box around
   xyres@lgLabelFontHeightF     = .02                ; label font height
   xyres@lgOrientation          = "vertical"

   xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
   xyres@vpWidthF    = 0.68     ; make the plot wider than it is high
   xyres@vpHeightF   = 0.4
   xyres@vpXF        = 0.11

   xyres@tiYAxisOffsetXF        = 0.015


   xyres@tiMainString="MOD, bbc2, pompa"


   xyres@xyDashPattern          = 0
   xyres@xyLineColor            = "black"
   xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
   xyres@lgLabelPosition        = "Left"
   xyres@pmLegendParallelPosF   = .8                  ; move units right
   xyres@pmLegendWidthF         = 0.08                ; Change width and
   xyres@pmLegendHeightF        = 0.15                ; height of legend.
   xyres@lgPerimOn              = False               ; turn off box around
   xyres@lgLabelFontHeightF     = .02                ; label font height
   xyres@lgOrientation          = "vertical"
   xyres@trYReverse             = False
   xyres@tiYAxisString          = "[mm h~S~-1~N~]"

   xyres@xyDashPattern          = 0
   xyres@trYMaxF                =  1.
   xyres@trYMinF                = -1.
   xyres@xyExplicitLegendLabels = "rain"
   xyres@pmLegendOrthogonalPosF = -1.35
   xyres@xyLineColor            = "blue"
   plot0  = gsn_csm_xy(wks,timeo,rm,xyres)
   ;xyres@xyDashPattern          = 1
   ;plot0   = gsn_csm_xy(wks,timeo,pm,res)
   ;overlay(plot0,plot)
   xyres@xyDashPattern          = 0
   xyres@xyExplicitLegendLabels = "evaporation"
   xyres@pmLegendOrthogonalPosF = -1.3
   xyres@xyLineColor            = "magenta"
   plot  = gsn_csm_xy(wks,timeo,embs+empl_sums,xyres)
   ;plot  = gsn_csm_xy(wks,timeo,ema,res)
   overlay(plot0,plot)
   xyres@xyExplicitLegendLabels = "surface runoff"
   xyres@pmLegendOrthogonalPosF = -1.25
   xyres@xyLineColor            = "green"
   plot  = gsn_csm_xy(wks,timeo,rsm,xyres)
   overlay(plot0,plot)
   xyres@xyExplicitLegendLabels = "ground runoff"
   xyres@pmLegendOrthogonalPosF = -1.2
   xyres@xyLineColor            = "orange"
   plot  = gsn_csm_xy(wks,timeo,rgm,xyres)
   overlay(plot0,plot)
   xyres@xyExplicitLegendLabels = "snow melt"
   xyres@pmLegendOrthogonalPosF = -1.15
   xyres@xyLineColor            = "cyan"
   plot  = gsn_csm_xy(wks,timeo,mm,xyres)
   overlay(plot0,plot)
   xyres@xyDashPattern          = 1
   xyres@xyExplicitLegendLabels = "sum"
   xyres@pmLegendOrthogonalPosF = -1.1
   xyres@xyLineColor            = "red"
   plot  = gsn_csm_xy(wks,timeo,sm,xyres)
   overlay(plot0,plot)
   xyres@xyDashPattern          = 3
   xyres@xyExplicitLegendLabels = "dW~B~SO~N~/dt"
   xyres@pmLegendOrthogonalPosF = -1.05
   xyres@xyLineColor            = "black"
   ;plot  = gsn_csm_xy(wks,timeo,dwdt_sum+dwidt+dwsdt,res)
   plot  = gsn_csm_xy(wks,timeo,dw,xyres)
   overlay(plot0,plot)
   res@xyDashPattern          = 2
   ;res@xyExplicitLegendLabels = "dW~B~SNOW~N~/dt"
   ;res@pmLegendOrthogonalPosF = -1.0
   ;res@xyLineColor            = "black"
   ;plot  = gsn_csm_xy(wks,timeo,dwsdt,res)
   ;overlay(plot0,plot)
   xyres@xyDashPattern          = 1
   xyres@xyExplicitLegendLabels = "dW~B~I~N~/dt"
   xyres@pmLegendOrthogonalPosF = -0.95
   xyres@xyLineColor            = "black"
   plot  = gsn_csm_xy(wks,timeo,dwi,xyres)
   overlay(plot0,plot)
   draw(plot0)
   frame(wks)
   print(dwi)
   xyres@xyDashPattern          = 0
   xyres@xyLineColor            = "black"
   xyres@pmLegendDisplayMode    = "Always"            ; turn on legend
   xyres@lgLabelPosition        = "Left"
   xyres@pmLegendParallelPosF   = .8                  ; move units right
   xyres@pmLegendWidthF         = 0.08                ; Change width and
   xyres@pmLegendHeightF        = 0.15                ; height of legend.
   xyres@lgPerimOn              = False               ; turn off box around
   xyres@lgLabelFontHeightF     = .02                ; label font height
   xyres@lgOrientation          = "vertical"
   xyres@trYReverse             = False
   xyres@tiYAxisString          = "[mm h~S~-1~N~]"

   xyres@xyDashPattern          = 0
   xyres@trYMaxF                = 10.
   xyres@trYMinF                = -5.
   xyres@xyExplicitLegendLabels = "T2m"
   xyres@pmLegendOrthogonalPosF = -1.35
   xyres@xyLineColor            = "blue"
   plot0  = gsn_csm_xy(wks,timeo,t2mm,xyres)
   draw(plot0)
   frame(wks)

   print(figname+".ps")

end
