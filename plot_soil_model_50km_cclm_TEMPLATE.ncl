load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; order of loading is important
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------------
begin

   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+SIMULATION+"/postprocessed/monmean/"
;   folder = "/project/pr04/regenass/cosmo5_validation/"+SIMULATION+"/monmean/"

   in     = addfile(folder+"monmean_out01.nc","r")
   s      = in->S_SO ; saturation
   wso      = in->W_SO ; saturation
   wio      = in->W_SO_ICE ; soil ice
   ;klw    = in->KLW ; conductivity at layer k-0.5
   ;dlw    = in->DLW ; conductivity at layer k-0.5
   zh     = in->soil1
   nz     = dimsizes(zh)
   level  = ispan(0,nz-1,1)
   ;in     = addfile(folder+"soil_michi.nc","r")
   sl     = in->SATLEV ; lowest level that is completely saturated
   wt     = in->WTDEPTH ; lowest level that is completely saturated
   ;q      = in->Q  ; flux

;   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+SIMULATION+"/postprocessed/monmean/"
   in     = addfile(folder+"monmean_out04.nc","r")
   evap_a  = in->ALHFL_S
   evap_bs = in->ALHFL_BS
   evap_pl = in->ALHFL_PL

;   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+SIMULATION+"/postprocessed/monmean/"
   in     = addfile(folder+"monmean_out05.nc","r")
   rs     = in->RUNOFF_S
   rg     = in->RUNOFF_G
   melt   = in->SNOW_MELT
   ;in     = addfile(folder+"soil.nc","r")
   ;melt   = in->SNOW_MELT

   time   = in->time
   timeo  = time/86400.; convert to days
   timeo  = timeo-365.-366.-15.5; convert to days
   nt     = dimsizes(time)
   print(timeo)
;   folder ="/scratch/snx3000/regenass/cosmo5_validation/"+SIMULATION+"/postprocessed/monmean/"
   in     = addfile(folder+"monmean_out03.nc","r")
   rain  = in->RAIN_CON+in->RAIN_GSP 
   ;prec  = in->TOT_PREC 
   ;raini  = in->RAIN_CON+in->RAIN_GSP 
   ;dewi   = in->DEW 
   ;preci  = in->TOT_PREC 
   timer  = in->time
   ntr    = dimsizes(timer)

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
   rain = rain * 6. ; convert to mm/6h
   do tt=0,nt-1
      rain(tt,:,:)   = where(ismissing(masks),-999.,rain(tt,:,:))
      ;dew(tt,:,:)   = where(ismissing(masks),-999.,dew(tt,:,:))
      ;prec(tt,:,:)   = where(ismissing(masks),-999.,prec(tt,:,:))
   end do
   rain@_FillValue = -999.
   ;dew@_FillValue = -999.
   ;prec@_FillValue = -999.

   do tt=0,nt-1
      evap_a(tt,:,:)   = where(ismissing(masks),-999.,evap_a(tt,:,:))
      evap_bs(tt,:,:)   = where(ismissing(masks),-999.,evap_bs(tt,:,:))
      do kk=0,nz-2
         evap_pl(tt,kk,:,:)   = where(ismissing(masks),-999.,evap_pl(tt,kk,:,:))
         wso(tt,kk,:,:)       = where(ismissing(masks),-999.,wso(tt,kk,:,:))
         wio(tt,kk,:,:)       = where(ismissing(masks),-999.,wio(tt,kk,:,:))
      end do
   end do
   evap_a@_FillValue = -999.
   evap_bs@_FillValue = -999.
   evap_pl@_FillValue = -999.
   wso@_FillValue = -999.
   wio@_FillValue = -999.

   sm   = dim_avg_n(s,(/2,3/))
   wm    = dim_avg_n(wso,(/2,3/))
   wim   = dim_avg_n(wio,(/2,3/))
   rm    = dim_avg_n(rain,(/1,2/))
   ;dm    = dim_avg_n(dew,(/1,2/))
   ;pm    = dim_avg_n(prec,(/1,2/))
   ema   = dim_avg_n(evap_a,(/1,2/))/2.51e06*21600.
   embs  = dim_avg_n(evap_bs,(/1,2/))/2.51e06*21600.
   empl  = dim_avg_n(evap_pl,(/2,3/))/2.51e06*21600.
   rsm   = dim_avg_n(rs,(/1,2/))
   rgm   = dim_avg_n(rg,(/1,2/))
   mm    = dim_avg_n(melt,(/1,2/))
   slm  = dim_avg_n(sl,(/1,2/))
   wtm  = dim_avg_n(wt,(/1,2/))

   rgm = rgm

   dwdt   = new((/nt,nz/),"float",-999)

   do tt=1,nt-1
      ;dwdt(tt) = dim_sum_n((w(tt,:)-w(tt-1,:))/(time(tt)-time(tt-1)),0)
      dwdt(tt,:) = (wm(tt,:)-wm(tt-1,:))/doubletofloat(time(tt)-time(tt-1))*1000.*3600.*6. ; convert to mm/6h
   end do

   ;dwdt_sum  = dim_sum(dwdt(:,0:7))
   ;w_sum     = dim_sum(wm(:,0:7))
   ;empl_sum  = dim_sum(empl(:,0:7))
   dwdt_sum  = dim_sum(dwdt(:,0:8))
   w_sum     = dim_sum(wm(:,0:8))
   wi_sum    = dim_sum(wim(:,0:8))
   empl_sum  = dim_sum(empl(:,0:8))
   sum_all = rm+rsm+rgm+ema+mm
   ;sum_all = rm+rsm+rgm+embs+empl_sum+mm

   rel_error = abs(dwdt_sum-sum_all)/rm*100.
   error = dwdt_sum-sum_all
   accum_error = error
   accum_error = 0.
   do tt=1,nt-1
      accum_error(tt) = accum_error(tt-1)+error(tt)
   end do
   accum_error = accum_error*4.*30.

   ; define resources

   res                        = True
   res@gsnFrame               = False ; don't advance frame
   res@gsnDraw                = False ; don't draw
   res@gsnMaximize            = True
   res@gsnPaperOrientation    = "landscape"

   res@tiXAxisString          = "Time [Days since 1980-01-01]"
   res@trYReverse             = True
   res@tmXBLabelsOn           = True
   res@tmXTLabelsOn           = False
   res@tmYRLabelsOn           = False
   res@tmYLLabelsOn           = True
   res@xyLineThicknessF       = 2
   res@pmLegendDisplayMode    = "Always"            ; turn on legend
   res@lgLabelPosition        = "Left"
   res@pmLegendParallelPosF   = .8                  ; move units right
   res@pmLegendWidthF         = 0.08                ; Change width and
   res@pmLegendHeightF        = 0.15                ; height of legend.
   res@lgPerimOn              = False               ; turn off box around
   res@lgLabelFontHeightF     = .02                ; label font height
   res@lgOrientation          = "vertical"

   res@pmLegendDisplayMode    = "Always"            ; turn on legend
   res@vpWidthF    = 0.68     ; make the plot wider than it is high
   res@vpHeightF   = 0.4
   res@vpXF        = 0.11

   figname = folder+"Budget_Plots/soil_model_"+SIMULATION
   wks     = gsn_open_wks("ps",figname)

   res@trYMaxF                = 119.
   res@trYMinF                = 0.

;;===================================================================

   xyres                      = True
   xyres@gsnFrame             = False ; don't advance frame
   xyres@gsnDraw              = False ; don't draw
   xyres@gsnMaximize            = True
   xyres@gsnPaperOrientation    = "landscape"

   xyres@tiXAxisString         = "Time [Years]"
   xyres@tmXBLabelsOn          = True
   xyres@tmXTLabelsOn          = False
   xyres@tmYRLabelsOn          = False
   xyres@tmYLLabelsOn          = True
   xyres@tmXBMode              = "Explicit"
   xyres@tmXBValues            = timeo(0::12)
   xyres@tmXBValues            = timeo(0:nt-1:12)
   xyres@tmXBMode              = "Explicit"
   xyres@tmXBLabels             = (/"1980-01","1981-01","1982-01","1983-01","1984-01","1985-01"/)
   xyres@tmXBLabelAngleF       = -45    ; tilt the XB labels 45 degrees

   xyres@cnLinesOn             = False
   xyres@cnFillOn              = True
   xyres@lbOrientation         = "vertical"
   xyres@trYReverse             = True

   xyres@cnFillColors          = (/-1,220,200,180,160,120,100,80,60,40,20/)
   xyres@cnLevelSelectionMode  = "ExplicitLevels"
   xyres@cnLevels              = fspan(0,100.,11)
   ;xyres@lbLabelStrings        =  (/"0.1","0.2","0.5","0.7","1.0","2.0","5.0","7.5","10","20"/)
   xyres@cnInfoLabelOn         = False
   xyres@cnExplicitLineLabelsOn= False

   xyres@vpWidthF    = 0.68     ; make the plot wider than it is high
   xyres@vpHeightF   = 0.35
   xyres@vpXF        = 0.11

   xyres@tiMainString=""
   xyres@tiMainOn=False
   res@pmLegendDisplayMode    = "Never"            ; turn on legend

   level@long_name = "soil layer"

   sm!0 = "Time"
   sm!1 = "level"
   sm&level=level
   sm&Time=timeo
   xyres@gsnRightString = "S [%]"

   plot  = gsn_csm_contour(wks,sm(level|:,Time|:),xyres)
   draw(plot)
   frame(wks)

   res@trXMaxF                = max(timeo)
   res@trXMinF                = 0.
   xyres@tmXBValues            = timeo(0::12)
   xyres@tmXBValues            = timeo(0:nt-1:12)
   xyres@tmXBMode              = "Explicit"
   xyres@tmXBLabels             = (/"1980-01","1981-01","1982-01","1983-01","1984-01","1985-01"/)
   res@tmXBLabelAngleF        = -45.    ; tilt the XB labels 45 degrees
   res@tmXBLabelJust          = "CenterLeft"
   res@xyLineColor            = "black"
   res@xyDashPattern          = 0
   res@trYMinF                = 1.
   res@trYMaxF                = 11.
   res@tiYAxisString          = "first fully saturated layer"
   plot  = gsn_csm_xy(wks,timeo,slm,res)
   draw(plot)
   frame(wks)

   res@trYMaxF                = 10.
   res@trYMinF                = -1.
   res@tiYAxisString          = "water table depth [m]"
   plot  = gsn_csm_xy(wks,timeo,wtm,res)
   draw(plot)
   frame(wks)

   res@trYReverse             = False
   res@trYMaxF                = 10.
   res@trYMinF                = 1.4
   res@tiYAxisString          = "soil water content [m]"
   res@xyLineColor            = "black"
   res@xyExplicitLegendLabels = "water"
   res@pmLegendOrthogonalPosF = -1.28
   plot  = gsn_csm_xy(wks,timeo,w_sum,res)
   res@xyExplicitLegendLabels = "ice"
   res@pmLegendOrthogonalPosF = -1.23
   res@xyDashPattern          = 2
   res@xyLineColor            = "orange"
   plot1 = gsn_csm_xy(wks,timeo,wi_sum,res)
   overlay(plot,plot1)
   draw(plot)
   frame(wks)

   res@xyDashPattern          = 0
   res@xyLineColor            = "black"
   res@pmLegendDisplayMode    = "Always"            ; turn on legend
   res@lgLabelPosition        = "Left"
   res@pmLegendParallelPosF   = 1.2                  ; move units right
   res@pmLegendWidthF         = 0.08                ; Change width and
   res@pmLegendHeightF        = 0.15                ; height of legend.
   res@lgPerimOn              = False               ; turn off box around
   res@lgLabelFontHeightF     = .02                ; label font height
   res@lgOrientation          = "vertical"
   res@trYReverse             = False
   res@tiYAxisString          = "[mm/6h]"

   res@xyDashPattern          = 0
   res@trYMaxF                =  1.
   res@trYMinF                = -1.
   res@xyExplicitLegendLabels = "rain"
   res@pmLegendOrthogonalPosF = -1.28
   res@xyLineColor            = "blue"
   plot0  = gsn_csm_xy(wks,timeo,rm,res)
   ;res@xyDashPattern          = 1
   ;plot   = gsn_csm_xy(wks,timeo,pm,res)
   ;overlay(plot0,plot)
   res@xyDashPattern          = 0
   res@xyExplicitLegendLabels = "evaporation"
   res@pmLegendOrthogonalPosF = -1.23
   res@xyLineColor            = "magenta"
   plot  = gsn_csm_xy(wks,timeo,ema,res)
   overlay(plot0,plot)
   res@xyExplicitLegendLabels = "surface runoff"
   res@pmLegendOrthogonalPosF = -1.18
   res@xyLineColor            = "green"
   plot  = gsn_csm_xy(wks,timeo,rsm,res)
   overlay(plot0,plot)
   res@xyExplicitLegendLabels = "ground runoff"
   res@pmLegendOrthogonalPosF = -1.13
   res@xyLineColor            = "orange"
   plot  = gsn_csm_xy(wks,timeo,rgm,res)
   overlay(plot0,plot)
   res@xyExplicitLegendLabels = "snow melt"
   res@pmLegendOrthogonalPosF = -1.08
   res@xyLineColor            = "cyan"
   plot  = gsn_csm_xy(wks,timeo,mm,res)
   overlay(plot0,plot)
   res@xyDashPattern          = 1
   res@xyExplicitLegendLabels = "sum"
   res@pmLegendOrthogonalPosF = -1.03
   res@xyLineColor            = "red"
   plot  = gsn_csm_xy(wks,timeo,sum_all,res)
   overlay(plot0,plot)
   res@xyDashPattern          = 3
   res@xyExplicitLegendLabels = "dW~B~SO~N~/dt"
   res@pmLegendOrthogonalPosF = -0.98
   res@xyLineColor            = "black"
   plot  = gsn_csm_xy(wks,timeo,dwdt_sum,res)
   overlay(plot0,plot)
   draw(plot0)
   frame(wks)

   res@pmLegendParallelPosF   = 1.2                  ; move units right
   res@xyDashPattern          = 0
   res@trYMaxF                =  0.8
   res@trYMinF                = -0.8
   res@xyExplicitLegendLabels = "rain"
   res@pmLegendOrthogonalPosF = -1.3
   res@xyLineColor            = "blue"
   plot0  = gsn_csm_xy(wks,timeo,rm,res)
   ;res@xyDashPattern          = 1
   ;res@pmLegendOrthogonalPosF = -1.25
   ;res@xyExplicitLegendLabels = "dew"
   ;res@xyLineColor            = "darkgreen"
   ;plot   = gsn_csm_xy(wks,timeo,dm_a,res)
   ;overlay(plot0,plot)
   res@xyDashPattern          = 0
   res@xyExplicitLegendLabels = "evap atm"
   res@pmLegendOrthogonalPosF = -1.2
   res@xyLineColor            = "grey"
   plot  = gsn_csm_xy(wks,timeo,ema,res)
   overlay(plot0,plot)
   res@pmLegendOrthogonalPosF = -1.15
   res@xyLineColor            = "magenta"
   res@xyExplicitLegendLabels = "evap soil"
   plot  = gsn_csm_xy(wks,timeo,embs+empl_sum,res)
   overlay(plot0,plot)
   res@xyDashPattern          = 1
   res@pmLegendOrthogonalPosF = -1.1
   res@xyExplicitLegendLabels = "bare soil evap"
   plot  = gsn_csm_xy(wks,timeo,embs,res)
   overlay(plot0,plot)
   res@pmLegendOrthogonalPosF = -1.05
   res@xyExplicitLegendLabels = "transpiration"
   res@xyDashPattern          = 2
   plot  = gsn_csm_xy(wks,timeo,empl_sum,res)
   overlay(plot0,plot)
   res@xyDashPattern          = 0
   res@xyExplicitLegendLabels = "surface runoff"
   res@pmLegendOrthogonalPosF = -1.0
   res@xyLineColor            = "green"
   plot  = gsn_csm_xy(wks,timeo,rsm,res)
   overlay(plot0,plot)
   res@xyExplicitLegendLabels = "ground runoff"
   res@pmLegendOrthogonalPosF = -0.95
   res@xyLineColor            = "orange"
   plot  = gsn_csm_xy(wks,timeo,rgm,res)
   overlay(plot0,plot)
   res@xyExplicitLegendLabels = "snow melt"
   res@pmLegendOrthogonalPosF = -0.9
   res@xyLineColor            = "cyan"
   plot  = gsn_csm_xy(wks,timeo,mm,res)
   overlay(plot0,plot)
   res@xyDashPattern          = 1
   res@xyExplicitLegendLabels = "sum"
   res@pmLegendOrthogonalPosF = -0.85
   res@xyLineColor            = "red"
   plot  = gsn_csm_xy(wks,timeo,sum_all,res)
   overlay(plot0,plot)
   res@xyDashPattern          = 3
   res@xyExplicitLegendLabels = "dW~B~SO~N~/dt"
   res@pmLegendOrthogonalPosF = -0.8
   res@xyLineColor            = "black"
   plot  = gsn_csm_xy(wks,timeo,dwdt_sum,res)
   overlay(plot0,plot)
   draw(plot0)
   frame(wks)

   res@pmLegendParallelPosF   = .8                  ; move units right
   res@xyDashPattern          = 0

   res@trYReverse             = False
   res@trYMaxF                = 100.
   res@trYMinF                = 0.
   res@tiYAxisString          = "relative error [%]"
   plot  = gsn_csm_xy(wks,timeo,rel_error,res)
   draw(plot)
   frame(wks)

   res@trYReverse             = False
   res@trYMaxF                = 100.
   res@trYMinF                = -100.
   res@tiYAxisString          = "accumulated error [mm]"
   plot  = gsn_csm_xy(wks,timeo,accum_error,res)
   draw(plot)
   frame(wks)


   ; PRUDENCE regions
   AL = (/44.,48.,5.,15./)    ; AL: Alps
   MD = (/36.,44.,3.,25./)    ; MD: Mediterranean
   BI = (/50.,59.,-10.,2./)   ; BI: British Isles
   IP = (/36.,44.,-10.,3./)   ; IP: Iberian Peninsula
   FR = (/44.,50.,-5.,5./)    ; FR: France
   ME = (/48.,55.,-2.,16./)   ; ME: Mid-Europe
   SC = (/55.,70.,5.,30./)    ; SC: Scandinavia
   EA = (/44.,55.,16.,30/)    ; EA: Eastern Europe
   ; additional regions of interest
   NO = (/69.,70.,25.,26/)    ; NO: Northern Norway


   rg_str = (/"AL","MD","BI","IP","FR","ME","SC","EA","NO"/)

   regions = (/AL,MD,BI,IP,FR,ME,SC,EA,NO/)
   S = dimsizes(regions)
   nr = S(0)

   ; loop over prudence regions
   do rr=0,nr-1
         
      RG = regions(rr,:)

      sr       = mask(s,(lat2d.ge.RG(0) .and. lat2d.le.RG(1).and. lon2d.ge.RG(2).and.lon2d.le.RG(3)), True)
      wsor     = mask(wso,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      wior     = mask(wio,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      rainr    = mask(rain,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      evap_bsr = mask(evap_bs,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      evap_plr = mask(evap_pl,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      rsr      = mask(rs,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      rgr      = mask(rg,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      meltr    = mask(melt,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      wtr      = mask(wt,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)
      slr      = mask(sl,(lat2d.ge.RG(0).and.lat2d.le.RG(1).and.lon2d.ge.RG(2).and.lon2d.le.RG(3)),True)

      sm    = dim_avg_n(sr,(/2,3/))
      wm    = dim_avg_n(wsor,(/2,3/))
      wim   = dim_avg_n(wior,(/2,3/))
      rm    = dim_avg_n(rainr,(/1,2/))
      embs  = dim_avg_n(evap_bsr,(/1,2/))/2.51e06*21600.
      empl  = dim_avg_n(evap_plr,(/2,3/))/2.51e06*21600.
      rsm   = -dim_avg_n(rsr,(/1,2/))
      rgm   = -dim_avg_n(rgr,(/1,2/))
      mm    = dim_avg_n(meltr,(/1,2/))
      slm   = dim_avg_n(slr,(/1,2/))
      wtm   = dim_avg_n(wtr,(/1,2/))


      dwdt   = new((/nt,nz/),"float",-999)

      do tt=1,nt-1
      ;dwdt(tt) = dim_sum_n((w(tt,:)-w(tt-1,:))/(time(tt)-time(tt-1)),0)
         dwdt(tt,:) = (wm(tt,:)-wm(tt-1,:))/doubletofloat(time(tt)-time(tt-1))*1000.*3600.*6. ; convert to mm/6h
      end do

      dwdt_sum  = dim_sum(dwdt(:,0:8))
      w_sum     = dim_sum(wm(:,0:8))
      wi_sum    = dim_sum(wim(:,0:8))
      empl_sum  = dim_sum(empl(:,0:8))
      sum_all = rm+rsm+rgm+embs+empl_sum+mm

      ;rel_error = abs(dwdt_sum-sum_all)/rm*100.
      error = dwdt_sum-sum_all
      accum_error = error
      accum_error = 0.
      do tt=1,nt-1
         accum_error(tt) = accum_error(tt-1)+error(tt)
      end do
      accum_error = accum_error*4.*30.


      figname = folder+"Budget_Plots/soil_model_"+SIMULATION+"_"+rg_str(rr)
      wks     = gsn_open_wks("ps",figname)


      level@long_name = "soil layer"

      sm!0 = "Time"
      sm!1 = "level"
      sm&level=level
      sm&Time=timeo
      xyres@pmLegendDisplayMode    = "Never"            ; turn on legend
      xyres@gsnRightString = "S [%]"
      xyres@gsnLeftString=rg_str(rr)

      plot  = gsn_csm_contour(wks,sm(level|:,Time|:),xyres)
      draw(plot)
      frame(wks)
      res@gsnMaximize  =  True
      res@gsnPaperOrientation = "landscape"
      res@vpXF = 1.5
      res@vpYF = 1.0
      res@pmLegendDisplayMode    = "Never"            ; turn on legend
      res@trYReverse             = True
      res@tiMainString=rg_str(rr)
      res@trXMaxF                = max(timeo);9855.
      res@trXMinF                = 0.
      xyres@tmXBValues            = timeo(0::12)
      xyres@tmXBValues            = timeo(0:nt-1:12)
      xyres@tmXBMode              = "Explicit"
      xyres@tmXBLabels             = (/"1980-01","1981-01","1982-01","1983-01","1984-01","1985-01"/)
      res@tmXBLabelAngleF        = -90    ; tilt the XB labels 45 degrees
      res@xyLineColor            = "black"
      res@xyDashPattern          = 0
      res@trYMinF                = 1.
      res@trYMaxF                = 10.
      res@tiYAxisString          = "first fully saturated layer"
      plot  = gsn_csm_xy(wks,timeo,slm,res)
      draw(plot)
      frame(wks)

      res@tiYAxisOffsetXF        = 0.015
      res@trYMaxF                = 6.
      res@trYMinF                = -1.
      res@tiYAxisString          = "water table depth [m]"
      plot  = gsn_csm_xy(wks,timeo,wtm,res)
      draw(plot)
      frame(wks)

      res@tiYAxisOffsetXF        = 0.0
      res@trYReverse             = False
      res@trYMaxF                = 10.0
      res@trYMinF                = 0.
      res@tiYAxisString          = "soil water content [m]"
      res@xyExplicitLegendLabels = "water"
      res@pmLegendOrthogonalPosF = -1.28
      plot  = gsn_csm_xy(wks,timeo,w_sum,res)
      res@xyDashPattern          = 2
      res@xyLineColor            = "orange"
      res@xyExplicitLegendLabels = "ice"
      res@pmLegendOrthogonalPosF = -1.23
      plot1 = gsn_csm_xy(wks,timeo,wi_sum,res)
      overlay(plot,plot1)
      draw(plot)
      frame(wks)

      res@xyDashPattern          = 0
      res@xyLineColor            = "black"
      res@pmLegendDisplayMode    = "Always"            ; turn on legend
      res@lgLabelPosition        = "Bottom"
      res@pmLegendParallelPosF   = .8                  ; move units right
      res@pmLegendWidthF         = 0.08                ; Change width and
      res@pmLegendHeightF        = 0.15                ; height of legend.
      res@lgPerimOn              = False               ; turn off box around
      res@lgLabelFontHeightF     = .02                ; label font height
      res@lgOrientation          = "vertical"
      res@trYReverse             = False
      res@tiYAxisString          = "[mm/6h]"

      res@xyDashPattern          = 0
      res@trYMaxF                =  2.
      res@trYMinF                = -2.
      res@xyExplicitLegendLabels = "rain"
      res@pmLegendOrthogonalPosF = -1.28
      res@xyLineColor            = "blue"
      plot0  = gsn_csm_xy(wks,timeo,rm,res)
      res@xyDashPattern          = 0
      res@xyExplicitLegendLabels = "evaporation"
      res@pmLegendOrthogonalPosF = -1.23
      res@xyLineColor            = "magenta"
      plot  = gsn_csm_xy(wks,timeo,ema,res)
      overlay(plot0,plot)
      res@xyExplicitLegendLabels = "surface runoff"
      res@pmLegendOrthogonalPosF = -1.18
      res@xyLineColor            = "green"
      plot  = gsn_csm_xy(wks,timeo,rsm,res)
      overlay(plot0,plot)
      res@xyExplicitLegendLabels = "ground runoff"
      res@pmLegendOrthogonalPosF = -1.13
      res@xyLineColor            = "orange"
      plot  = gsn_csm_xy(wks,timeo,rgm,res)
      overlay(plot0,plot)
      res@xyExplicitLegendLabels = "snow melt"
      res@pmLegendOrthogonalPosF = -1.08
      res@xyLineColor            = "cyan"
      plot  = gsn_csm_xy(wks,timeo,mm,res)
      overlay(plot0,plot)
      res@xyDashPattern          = 1
      res@xyExplicitLegendLabels = "sum"
      res@pmLegendOrthogonalPosF = -1.03
      res@xyLineColor            = "red"
      plot  = gsn_csm_xy(wks,timeo,sum_all,res)
      overlay(plot0,plot)
      res@xyDashPattern          = 3
      res@xyExplicitLegendLabels = "dW~B~SO~N~/dt"
      res@pmLegendOrthogonalPosF = -0.98
      res@xyLineColor            = "black"
      plot  = gsn_csm_xy(wks,timeo,dwdt_sum,res)
      overlay(plot0,plot)
      draw(plot0)
      frame(wks)

      res@trYReverse             = False
      res@trYMaxF                = 100.
      res@trYMinF                = -100.
      res@tiYAxisString          = "accumulated error [mm]"
      plot  = gsn_csm_xy(wks,timeo,accum_error,res)
      draw(plot)
      frame(wks)


   end do

   print(figname+".ps")

end
