<script>
  
  import { scaleLinear } from "d3-scale";
  // import { Date } from "d3-time"
  import Chart from './Chart.svelte';
  import { onMount } from 'svelte';
  import { selectAll } from 'd3-selection'
  import { drag } from 'd3-drag';
  import queryString from "query-string";
  import Checkbox from './Checkbox.svelte';
  import Arrow from './Arrow.svelte';
  import { format } from 'd3-format'
  import { event } from 'd3-selection'

  import katex from 'katex';

  const legendheight = 67 

  const legendheight_hosp = 40 

  function range(n){
    return Array(n).fill().map((_, i) => i);
  }

  function formatNumber(num) {
    return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
  }

  var sum = function(arr, bools){
    var x = 0
    for (var i = 0; i < arr.length; i++) {
      x = x + arr[i]*(bools[i] ? 1 : 0)
    }
    return x
  }

  var Integrators = {
    Euler    : [[1]],
    Midpoint : [[.5,.5],[0, 1]],
    Heun     : [[1, 1],[.5,.5]],
    Ralston  : [[2/3,2/3],[.25,.75]],
    K3       : [[.5,.5],[1,-1,2],[1/6,2/3,1/6]],
    SSP33    : [[1,1],[.5,.25,.25],[1/6,1/6,2/3]],
    SSP43    : [[.5,.5],[1,.5,.5],[.5,1/6,1/6,1/6],[1/6,1/6,1/6,1/2]],
    RK4      : [[.5,.5],[.5,0,.5],[1,0,0,1],[1/6,1/3,1/3,1/6]],
    RK38     : [[1/3,1/3],[2/3,-1/3,1],[1,1,-1,1],[1/8,3/8,3/8,1/8]]
  };

  // f is a func of time t and state y
  // y is the initial state, t is the time, h is the timestep
  // updated y is returned.
  var integrate=(m,f,y,t,h)=>{
    for (var k=[],ki=0; ki<m.length; ki++) {
      var _y=y.slice(), dt=ki?((m[ki-1][0])*h):0;
      for (var l=0; l<_y.length; l++) for (var j=1; j<=ki; j++) _y[l]=_y[l]+h*(m[ki-1][j])*(k[ki-1][l]);
      k[ki]=f(t+dt,_y,dt); 
    }
    for (var r=y.slice(),l=0; l<_y.length; l++) for (var j=0; j<k.length; j++) r[l]=r[l]+h*(k[j][l])*(m[ki-1][j]);
    return r;
  }


  $: Time_to_death     = 21
  $: logN              = Math.log(10500000)
  $: N                 = Math.exp(logN)
  $: I0                = 1
  $: R0                = 1.4
  $: D_incbation       = 5.2       
  $: D_infectious      = 2.9 
  $: D_recovery_mild   = (8 - 2.9)  
  $: D_recovery_severe = 6.9
  $: D_hospital_lag    = 5
  $: D_death           = Time_to_death - D_infectious 
  $: CFR               = 0.008  
  $: InterventionTime  = 22
  $: OMInterventionAmt = -0.35
  $: InterventionAmt   = 1 + OMInterventionAmt
  $: Time              = 220
  $: Xmax              = 110000
  $: dt                = 2
  $: P_SEVERE          = 0.04
  $: duration          = 7*12*1e10
  $: startSimulationAtDeaths = 120

  $: state = location.protocol + '//' + location.host + location.pathname + "?" + queryString.stringify({"startSimulationAtDeaths":startSimulationAtDeaths,              
               "R0":R0,              
               "InterventionAmt":InterventionAmt
              })

// dt, N, I0, R0, D_incbation, D_infectious, D_recovery_mild, D_hospital_lag, D_recovery_severe, D_death, P_SEVERE, CFR, InterventionTime, InterventionAmt, duration

  function get_solution(dt, N, I0, R0, D_incbation, D_infectious, D_recovery_mild, D_hospital_lag, D_recovery_severe, D_death, P_SEVERE, CFR, InterventionTime, InterventionAmt, duration, startSimulationAtDeaths) {

    var interpolation_steps = 40
    var steps = 100*interpolation_steps
    var dt = dt/interpolation_steps
    var sample_step = interpolation_steps

    var method = Integrators["RK4"]

    

    InterventionTime = 1000

    function f(t, x){

      // SEIR ODE
      if (t > InterventionTime && t < InterventionTime + duration){
        var beta = (InterventionAmt)*R0/(D_infectious)
      } else if (t > InterventionTime + duration) {
        var beta = 0.5*R0/(D_infectious)        
      } else {
        var beta = R0/(D_infectious)
      }
      var a     = 1/D_incbation
      var gamma = 1/D_infectious
      
      var S        = x[0] // Susectable
      var E        = x[1] // Exposed
      var I        = x[2] // Infectious 
      var Mild     = x[3] // Recovering (Mild)     
      var Severe   = x[4] // Recovering (Severe at home)
      var Severe_H = x[5] // Recovering (Severe in hospital)
      var Fatal    = x[6] // Recovering (Fatal)
      var R_Mild   = x[7] // Recovered
      var R_Severe = x[8] // Recovered
      var R_Fatal  = x[9] // Dead

      var p_severe = P_SEVERE
      var p_fatal  = CFR
      var p_mild   = 1 - P_SEVERE - CFR

      var dS        = -beta*I*S
      var dE        =  beta*I*S - a*E
      var dI        =  a*E - gamma*I
      var dMild     =  p_mild*gamma*I   - (1/D_recovery_mild)*Mild
      var dSevere   =  p_severe*gamma*I - (1/D_hospital_lag)*Severe
      var dSevere_H =  (1/D_hospital_lag)*Severe - (1/D_recovery_severe)*Severe_H
      var dFatal    =  p_fatal*gamma*I  - (1/D_death)*Fatal
      var dR_Mild   =  (1/D_recovery_mild)*Mild
      var dR_Severe =  (1/D_recovery_severe)*Severe_H
      var dR_Fatal  =  (1/D_death)*Fatal
      var dI_abs    =  a*E

      //      0   1   2   3      4        5          6       7        8          9         10
      return [dS, dE, dI, dMild, dSevere, dSevere_H, dFatal, dR_Mild, dR_Severe, dR_Fatal, dI_abs]
    }


    
    var initial_v = [1 - I0/N, 0, I0/N, 0, 0, 0, 0, 0, 0, 0, 0]
    var initial_t = 0
    var start_t = 0

    var initial_steps = steps * 10
    while (initial_steps--) { 
      if ((initial_steps+1) % (sample_step) == 0) {
        var d_v = f(initial_t, initial_v)            
        if (N*d_v[9] > startSimulationAtDeaths)
        {       
          start_t = initial_t - (D_death + D_infectious) - D_incbation
          InterventionTime = initial_t
          //console.log('start_t', start_t, (N*d_v[9]))
          //console.log('InterventionTime', InterventionTime)
          break
        }
        
        // console.log((v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7] + v[8] + v[9]))
        // console.log(v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , v[6] , v[7] , v[8] , v[9])
      }
      initial_v =integrate(method,f,initial_v, initial_t,dt); 
      initial_t+=dt
    }


    var v = [1 - I0/N, 0, I0/N, 0, 0, 0, 0, 0, 0, 0, 0]
    var t = 0

    var P  = []
    var chart_P = []
    var chart_Hosp = []
    var TI = []
    var Iters = []
    steps += Math.floor(start_t / dt)

    //Age distribution of hospitalized patients in percents, taken from UZIS data 17.10.2020
    var age_distribution_hospitalized = [1.14, 2.7, 3.87, 8.85, 11.39, 18.12, 14.68, 14.55, 11.41, 13.28]

    var prevInfectious = 0
    while (steps--) { 
      if ((steps+1) % (sample_step) == 0) {
            //    Dead   Hospital          Recovered        Infectious   Exposed
        if (t > start_t)
        {
          P.push([ N*v[9], N*(v[5]+v[6]),  N*(v[7] + v[8]), N*v[2],    N*v[1] ])
          chart_P.push([ N*v[9], N*(v[5]+v[6]),  N*(v[7] + v[8]), v[10]*N - prevInfectious,    N*v[1] ])
          //console.log(v[10]*N)

          var hosp = N*(v[5]+v[6])
          var hosp_by_age = []
          for (var i=0; i < age_distribution_hospitalized.length; i++)
          {
            hosp_by_age.push(age_distribution_hospitalized[i] / 100.0 * hosp)    
          }
          //console.log(hosp_by_age)
          chart_Hosp.push(hosp_by_age)

                    
          Iters.push(v)
          TI.push(N*(1-v[0]))
          
        }
        
        prevInfectious =  v[10]*N
        // console.log((v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7] + v[8] + v[9]))
        // console.log(v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , v[6] , v[7] , v[8] , v[9])
      }
      v =integrate(method,f,v,t,dt); 
      t+=dt
    }


    return {"P": P, 
            "chart_P": chart_P, 
            "chart_Hosp": chart_Hosp, 
            "deaths": N*v[6], 
            "total": 1-v[0],
            "total_infected": TI,
            "Iters":Iters,
            "dIters": f}
  }

  function max(P, checked) {
    return P.reduce((max, b) => Math.max(max, sum(b, checked) ), sum(P[0], checked) )
  }

  $: Sol            = get_solution(dt, N, I0, R0, D_incbation, D_infectious, D_recovery_mild, D_hospital_lag, D_recovery_severe, D_death, P_SEVERE, CFR, InterventionTime, InterventionAmt, duration, startSimulationAtDeaths)
  $: P              = Sol["P"].slice(0,100)
  $: chart_P        = Sol["chart_P"].slice(0,100)
  $: chart_Hosp     = Sol["chart_Hosp"].slice(0,100)
  $: timestep       = dt
  $: tmax           = dt*100
  $: deaths         = Sol["deaths"]
  $: total          = Sol["total"]
  $: total_infected = Sol["total_infected"].slice(0,100)
  $: Iters          = Sol["Iters"]  
  $: dIters         = Sol["dIters"]
  $: Pmax           = max(chart_P, checked)
  $: Hospmax        = 40000
  $: lock           = false

  var colors = [ "#555555", "#8da0cb", "#4daf4a", "#f0027f", "#fdc086"]
  
  var colors_hosp = ["#f6fedb","#e6d3a3","#d8d174","#b6c454","#91972a","#7f7caf","#28587b","#644536","#b2675e","#ef6f6c"]

  var Plock = 1
  


  $: parsed = "";
  onMount(async () => {

    if (typeof window !== 'undefined') {
      parsed = queryString.parse(window.location.search)
      if (!(parsed.startSimulationAtDeaths === undefined)) {startSimulationAtDeaths = parsed.startSimulationAtDeaths}      
      if (!(parsed.R0 === undefined)) {R0 = parseFloat(parsed.R0)}      
      if (!(parsed.InterventionAmt === undefined)) {InterventionAmt = parseFloat(parsed.InterventionAmt)}   
    }
  });

  function lock_yaxis(){
    Plock = Pmax
    lock  = true
  }

  function unlock_yaxis(){
    lock = false
  }

  const padding = { top: 20, right: 0, bottom: 20, left: 25 };

  let width  = 750;
  let height = 400;

  $: xScaleTime = scaleLinear()
    .domain([0, tmax])
    .range([padding.left, width - padding.right]);

  $: xScaleTimeInv = scaleLinear()
    .domain([0, width])
    .range([0, tmax]);

  $: indexToTime = scaleLinear()
    .domain([0, P.length])
    .range([-25, tmax])

  window.addEventListener('mouseup', unlock_yaxis);

  $: current_day = 11
  $: checked = [true, true, false, true, false]
  $: active  = 0
  $: active_ = active >= 0 ? active : current_day

  $: checked_hosp = [true, true, true, true, true, true, true, true, true, true, true]
  $: active_hosp  = 0
  $: active_hosp_ = active_hosp >= 0 ? active_hosp : current_day

  var Tinc_s = "\\color{#CCC}{T^{-1}_{\\text{inc}}} "
  var Tinf_s = "\\color{#CCC}{T^{-1}_{\\text{inf}}}"
  var Rt_s   = "\\color{#CCC}{\\frac{\\mathcal{R}_{t}}{T_{\\text{inf}}}} "
  $: ode_eqn = katex.renderToString("\\frac{d S}{d t}=-" +Rt_s +"\\cdot IS,\\qquad \\frac{d E}{d t}=" +Rt_s +"\\cdot IS- " + Tinc_s + " E,\\qquad \\frac{d I}{d t}=" + Tinc_s + "E-" + Tinf_s+ "I, \\qquad \\frac{d R}{d t}=" + Tinf_s+ "I", {
    throwOnError: false,
    displayMode: true,
    colorIsTextColor: true
  });

  function math_inline(str) {
    return katex.renderToString(str, {
    throwOnError: false,
    displayMode: false,
    colorIsTextColor: true
    });
  }

  function math_display(str) {
    return katex.renderToString(str, {
    throwOnError: false,
    displayMode: true,
    colorIsTextColor: true
    });
  }
  
  $: p_num_ind = 40

  $: get_d = function(i){
    return dIters(indexToTime(i), Iters[i])
  }

  function get_milestones(P){

    function argmax(x, index) {
      return x.map((x, i) => [x[index], i]).reduce((r, a) => (a[0] > r[0] ? a : r))[1];
    }

     //    Dead   Hospital          Recovered        Infectious   Exposed
    var milestones = []    

    var i = argmax(P, 1)
    milestones.push([i*dt, "Maximum hospitalizací: " + format(",")(Math.round(P[i][1]))])

    return milestones
  }

  $: milestones = get_milestones(P)
  $: log = true

</script>

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.css" integrity="sha384-bsHo4/LA+lkZv61JspMDQB9QP1TtO4IgOf2yYS+J6VdAYLVyx1c3XKcsHh0Vy8Ws" crossorigin="anonymous">

<style>
  .small { font: italic 6px Source Code Pro; }
  @import url('https://fonts.googleapis.com/css?family=Source+Code+Pro&display=swap');


  h2 {
    margin: auto;
    width: 950px;
    font-size: 40px;
    padding-top: 20px;
    padding-bottom: 20px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    padding-bottom: 30px
  }

  .center {
    margin: auto;
    width: 950px;
    padding-bottom: 20px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    color:#666;
    font-size: 16.5px;
    text-align: justify;
    line-height: 24px
  }

  .ack {
    margin: auto;
    width: 950px;
    padding-bottom: 20px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    color:#333;
    font-size: 13px;
  }

  .row {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    margin: auto;
    display: flex;
    width: 948px;
    font-size: 13px;
  }

  .caption {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    font-size: 13px;    
  }

  .column {
    flex: 158px;
    padding: 0px 5px 5px 0px;
    margin: 0px 5px 5px 5px;
    /*border-top: 2px solid #999*/
  }

  .minorTitle {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    margin: auto;
    display: flex;
    width: 950px;
    font-size: 17px;
    color: #666;
  }

  .minorTitleColumn{
    flex: 60px;
    padding: 3px;
    border-bottom: 2px solid #999;
  }


  .paneltext{
    position:relative;
    height:130px;
  }

  .paneltitle{
    color:#777; 
    line-height: 17px; 
    padding-bottom: 4px;
    font-weight: 700;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
  }

  .paneldesc{
    color:#888; 
    text-align: left;
    font-weight: 300;
  }

  .slidertext{
    color:#555; 
    line-height: 7px; 
    padding-bottom: 0px; 
    padding-top: 7px;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    font-family: 'Source Code Pro', monospace;
    font-size: 10px;
    text-align: right;
    /*font-weight: bold*/
  }
    
  .range {
    width: 100%;
  }

  .chart {
    width: 100%;
    margin: 0 auto;
    padding-top:0px;
    padding-bottom:10px;
  }

  .legend {
    color: #888;
    font-family: Helvetica, Arial;
    font-size: .725em;
    font-weight: 200;
    height: 100px;
    left: 20px;
    top: 4px;
    position: absolute;
  }

  .legendtitle {
    color:#777; 
    font-size:13px;
    padding-bottom: 6px;
    font-weight: 600;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
  }


  .legendtext{
    color:#888; 
    font-size:13px;
    padding-bottom: 5px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    line-height: 14px;
  }

  .legendtextnum{
    color:#888; 
    font-size:13px;
    padding-bottom: 5px;
    font-weight: 300;
    line-height: 12px;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    left: -3px;
    position: relative;
  }

  .tick {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    font-size: .725em;
    font-weight: 200;
    font-size: 13px
  }

  td { 
    text-align: left;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    border-bottom: 1px solid #DDD;
    border-collapse: collapse;
    padding: 3px;
    /*font-size: 14px;*/
  }

  tr {
    border-collapse: collapse;
    border-spacing: 15px;
  }

  .eqn {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    margin: auto;
    display: flex;
    flex-flow: row wrap;
    width: 950px;
    column-count: 4;
    font-weight: 300;
    color:#666;
    font-size: 16.5px;
  }

  th { font-weight: 500; text-align: left; padding-bottom: 5px; vertical-align: text-top;     border-bottom: 1px solid #DDD; }

  a:link { color: grey; }
  a:visited { color: grey; }

</style>

<h2>Epidemická kalkulačka</h2>

<div class="minorTitle" >
  <div class="minorTitleColumn" style="margin-bottom:20px">Vývoj epidemie</div>        
</div>

<div class="chart" style="display: flex; max-width: 1120px">

  <div style="flex: 0 0 270px; width:270px;">
    <div style="position:relative; top:48px; right:-115px">
      <div class="legendtext" style="position:absolute; left:-16px; top:-34px; width:50px; height: 100px; font-size: 13px; line-height:16px; font-weight: normal; text-align: center"><b>Den</b><br> {Math.round(indexToTime(active_))}</div>

      <!-- Susceptible -->
      <div style="position:absolute; left:0px; top:0px; width: 180px; height: 100px">

        <span style="pointer-events: none"><Checkbox color="#CCC"/></span>
        <Arrow height="41"/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Náchylní infekci</div>
          <div style="padding-top: 5px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*Iters[active_][0]))} 
                                  ({ (100*Iters[active_][0]).toFixed(2) }%)</i></div>
          <div class="legendtextnum" style='visibility: {(active_ < (Iters.length - 1)) ? 'visible':'hidden'};'><span style="font-size:12px; padding-right:2px; color:#CCC;">Δ</span> <i>{formatNumber(Math.round(N*get_d(active_)[0]))} / den</i>
                                 </div>
          </div>
        </div>
          <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4px; position:relative;">Lidé, kteří se virem ještě nenakazili.</div>

      </div>

      <!-- Exposed -->
      <div style="position:absolute; left:0px; top:{legendheight*1}px; width: 180px; height: 100px">

        <Checkbox color="{colors[4]}" bind:checked={checked[4]}/>      
        <Arrow height="41"/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Nakažení</div>

          <div style="padding-top: 5px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*Iters[active_][1]))} 
                                  ({ (100*Iters[active_][1]).toFixed(2) }%)</div>
          <div class="legendtextnum" style='visibility: {(active_ < (Iters.length - 1)) ? 'visible':'hidden'};'><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*get_d(active_)[1])) } / den</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4px; position:relative;">Lidé v inkubační době.</div>

      </div>

      <!-- Infectious -->
      <div style="position:absolute; left:0px; top:{legendheight*2}px; width: 180px; height: 100px">

        <Checkbox color="{colors[3]}" bind:checked={checked[3]}/>
        <Arrow height="41"/>   

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Nově pozitivní</div>
          <div style="padding-top: 5px; padding-bottom: 1px">          
          <div class="legendtextnum" style='visibility: {(active_ < (Iters.length - 1)) ? 'visible':'hidden'};'><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(chart_P[active_][3])) } / den</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4px; position:relative;">Počet infekčních lidí vč. bezpříznakových.</div>


      </div>

      <!-- Removed -->
     <!-- <div style="position:absolute; left:0px; top:{legendheight*3}px; width: 180px; height: 100px">

        <Checkbox color="grey" callback={(s) => {checked[1] = s; checked[0] = s; checked[2] = s} }/>
        <Arrow height="56" arrowhead="" dasharray="3 2"/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Removed</div>
          <div style="padding-top: 10px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N* (1-Iters[active_][0]-Iters[active_][1]-Iters[active_][2]) ))} 
                                  ({ ((100*(1-Iters[active_][0]-Iters[active_][1]-Iters[active_][2]))).toFixed(2) }%)</div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[3]+get_d(active_)[4]+get_d(active_)[5]+get_d(active_)[6]+get_d(active_)[7]) )) } / den</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4x; position:relative;">Population no longer infectious due to isolation or immunity.</div>

      </div> -->

      <!-- Recovered -->
      <div style="position:absolute; left:0px; top:{legendheight*3}px; width: 180px; height: 100px">
        <Checkbox color="{colors[2]}" bind:checked={checked[2]}/>
        <!--<Arrow height="23" arrowhead="" dasharray="3 2"/>-->
        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Uzdravení</div>

          <div style="padding-top: 3px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][7]+Iters[active_][8]) ))} 
                                  ({ (100*(Iters[active_][7]+Iters[active_][8])).toFixed(2) }%)</div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 8px; position:relative;">Lidé imunní minimálně do skončení této vlny epidemie.</div>

      </div>

      <!-- Hospitalized -->
      <div style="position:absolute; left:0px; top:{legendheight*4}px; width: 180px; height: 100px">
        <Arrow height="43" arrowhead="" dasharray="3 2"/>
        <Checkbox color="{colors[1]}" bind:checked={checked[1]}/>
        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Hospitalizováno</div>
          <div style="padding-top: 3px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][5]+Iters[active_][6]) ))} 
                                  ({ (100*(Iters[active_][5]+Iters[active_][6])).toFixed(2) }%)</div>
          </div>
          <div class="legendtextnum" style='visibility: {(active_ < (Iters.length - 1)) ? 'visible':'hidden'};'><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[5]+get_d(active_)[6]))) } / den</i>
                                 </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 10px; position:relative;">Počet lidí aktuálně v nemocnicích.</div>

      </div>

      <div style="position:absolute; left:0px; top:{legendheight*5}px; width: 180px; height: 100px">
        <Arrow height="40" arrowhead="" dasharray="3 2"/>

        <Checkbox color="{colors[0]}" bind:checked={checked[0]}/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Zemřelí</div>
          <div style="padding-top: 3px; padding-bottom: 1px">          
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*Iters[active_][9]))} 
                                  ({ (100*Iters[active_][9]).toFixed(2) }%)</div>
          <div class="legendtextnum" style='visibility: {(active_ < (Iters.length - 1)) ? 'visible':'hidden'};'><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*get_d(active_)[9])) } / den</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; color:#ff0000; width:105px; left:-111px; top: 10px; position:relative">Zemřelí. Platí pouze za předpokladu, že se všem dostane lékařské péče (viz níže).</div>
      </div>
    </div>
  </div>

  <div style="flex: 0 0 890px; width:890px; height: {height+128}px; position:relative;">

    <div style="position:relative; top:60px; left: 10px">
      <Chart bind:checked={checked}
             bind:active={active}
             y = {chart_P} 
             xmax = {Xmax} 
             total_infected = {total_infected} 
             deaths = {deaths} 
             total = {total} 
             timestep={timestep}
             tmax={tmax}
             N={N}
             ymax={lock ? Plock: Pmax}
             InterventionTime={InterventionTime}
             colors={colors}
             log={!log}/>
      </div>

     
      <!-- Intervention Line -->
      <div style="position: absolute; width:{width+15}px; height: {height}px; position: absolute; top:100px; left:10px; pointer-events: none">
        <div id="dottedline"  style="pointer-events: all;
                    position: absolute;
                    top:-38px;
                    left:{xScaleTime(InterventionTime)}px;
                    visibility: {(xScaleTime(InterventionTime) < (width - padding.right)) ? 'visible':'hidden'};
                    width:2px;
                    background-color:#FFF;
                    border-right: 1px dashed black;
                    pointer-events: all;                    
                    height:{height+19}px">

        <div style="position:absolute; opacity: 0.5; top:-5px; left:10px; width: 120px">
        <span style="font-size: 13px">{@html math_inline("\\mathcal{R}=" + (R0*InterventionAmt).toFixed(2) )}</span> ⟶ 
        </div>

        {#if xScaleTime(InterventionTime) >= 100}
          <div style="position:absolute; opacity: 0.5; top:-2px; left:-97px; width: 120px">
          <span style="font-size: 13px">⟵ {@html math_inline("\\mathcal{R}=" + (R0).toFixed(2) )}</span>
          </div>      
        {/if}

        <div id="interventionDrag" class="legendtext" style="flex: 0 0 160px; width:120px; position:relative;  top:-70px; height: 60px; padding-right: 15px; left: -125px; pointer-events: all;" >
          <div class="paneltitle" style="top:9px; position: relative; text-align: right">Vývoj do dnešního dne</div>
          <span></span>          
        </div>


       


        </div>
      </div>

      <!-- Intervention Line slider -->
      <div style="position: absolute; width:{width+15}px; height: {height}px; position: absolute; top:120px; left:10px; pointer-events: none">
        <div style="
            position: absolute;
            top:-38px;
            left:{xScaleTime(InterventionTime)}px;
            visibility: {(xScaleTime(InterventionTime) < (width - padding.right)) ? 'visible':'hidden'};
            width:2px;
            background-color:#FFF;
            border-right: 1px dashed black;
            cursor:col-resize;
            height:{height}px">
            <div style="flex: 0 0 160px; width:300px; position:relative; top:-125px; left: 1px" >
              <div class="caption" style="pointer-events: none; position: absolute; left:0; top:40px; width:150px; border-left: 2px solid #777; padding: 5px 7px 7px 7px; ">      
              <div class="paneltext"  style="height:20px; text-align: right">
              <div class="paneldesc">snížit (-) nebo zvýšit (+) přenos viru o <br></div>
              </div>
              <div style="pointer-events: all">
              <div class="slidertext" on:mousedown={lock_yaxis}>{-(100*(1-InterventionAmt)).toFixed(2)}%</div>
              <input class="range" type=range bind:value={OMInterventionAmt} min=-0.5 max=0.5 step=0.01 on:mousedown={lock_yaxis}>
              </div>
              </div>
            </div>
          </div>
      </div>

      <div style="pointer-events: none;
                  position: absolute;
                  top:{height+84}px;
                  left:{0}px;
                  width:{780}px;
                  opacity: 1.0;
                  height:25px;
                  cursor:col-resize">
            {#each milestones as milestone}
              <div style="position:absolute; left: {xScaleTime(milestone[0])+8}px; top: -30px;">
                <span style="opacity: 0.3"><Arrow height=30 arrowhead="#circle" dasharray = "2 1"/></span>
                  <div class="tick" style="position: relative; left: 0px; top: 35px; max-width: 130px; color: #BBB; background-color: white; padding-left: 4px; padding-right: 4px">{@html milestone[1]}</div>
              </div>
            {/each}
      </div>       

   </div>

</div>

<div class="minorTitle">
  <div class="minorTitleColumn">Vývoj hospitalizací</div>        
</div>

<div class="chart" style="display: flex; max-width: 1120px">

  <div style="flex: 0 0 270px; width:270px;">
    <div style="position:relative; top:48px; right:-115px">
      <div class="legendtext" style="position:absolute; left:-16px; top:-34px; width:50px; height: 100px; font-size: 13px; line-height:16px; font-weight: normal; text-align: center"><b>Den</b><br> {Math.round(indexToTime(active_hosp_))}</div>
      
      <div style="position:absolute; left:0px; top:0px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[0]}" bind:checked={checked_hosp[0]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">0-17 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][0]))} ({ (100*chart_Hosp[active_hosp_][0] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*1}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[1]}" bind:checked={checked_hosp[1]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">18-29 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][1]))} ({ (100*chart_Hosp[active_hosp_][1] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*2}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[2]}" bind:checked={checked_hosp[2]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">30-39 let</div>
          <div style="padding-top: 1px; padding-bottom:5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][2]))} ({ (100*chart_Hosp[active_hosp_][2] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*3}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[3]}" bind:checked={checked_hosp[3]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">40-49 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][3]))} ({ (100*chart_Hosp[active_hosp_][3] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*4}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[4]}" bind:checked={checked_hosp[4]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">50-59 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][4]))} ({ (100*chart_Hosp[active_hosp_][4] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*5}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[5]}" bind:checked={checked_hosp[5]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">60-69 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][5]))} ({ (100*chart_Hosp[active_hosp_][5] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*6}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[6]}" bind:checked={checked_hosp[6]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">70-74 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][6]))} ({ (100*chart_Hosp[active_hosp_][6] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*7}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[7]}" bind:checked={checked_hosp[7]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">75-79 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][7]))} ({ (100*chart_Hosp[active_hosp_][7] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*8}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[8]}" bind:checked={checked_hosp[8]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">79-84 let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][8]))} ({ (100*chart_Hosp[active_hosp_][8] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>
      <div style="position:absolute; left:0px; top:{legendheight_hosp*9}px; width: 180px; height: 100px">
        <Checkbox color="{colors_hosp[9]}" bind:checked={checked_hosp[9]}/>
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">85+  let</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(chart_Hosp[active_hosp_][9]))} ({ (100*chart_Hosp[active_hosp_][9] / sum(chart_Hosp[active_hosp_], checked_hosp) ).toFixed(2) }%) </i></div>          
          </div>
        </div>          
      </div>      

      <div style="position:absolute; left:0px; top:{legendheight_hosp*10+10}px; width: 180px; height: 100px">   
        <div class="legend" style="position:absolute;">          
          <div class="legendtitle">Celkem</div>
          <div style="padding-top: 1px; padding-bottom: 5px">
          <div class="legendtextnum"><i>{formatNumber(Math.round(sum(chart_Hosp[active_hosp_], checked_hosp)))}  </i></div>          
          </div>
        </div>          
      </div>    
      
    </div>
  </div>
  <div style="flex: 0 0 890px; width:890px; height: {height+128}px; position:relative;">

    <div style="position:relative; top:60px; left: 10px">
      <Chart bind:checked={checked_hosp}
             bind:active={active_hosp}
             y = {chart_Hosp} 
             xmax = {Xmax} 
             total_infected = {total_infected} 
             deaths = {deaths} 
             total = {total} 
             timestep={timestep}
             tmax={tmax}
             N={N}
             ymax={lock ? Plock: Hospmax}
             InterventionTime={InterventionTime}
             colors={colors_hosp}
             log={!log}
             show_capacity={true}/>
      </div>

     
      <!-- Intervention Line -->
      <div style="position: absolute; width:{width+15}px; height: {height}px; position: absolute; top:100px; left:10px; pointer-events: none">
        <div id="dottedline"  style="pointer-events: all;
                    position: absolute;
                    top:-38px;
                    left:{xScaleTime(InterventionTime)}px;
                    visibility: {(xScaleTime(InterventionTime) < (width - padding.right)) ? 'visible':'hidden'};
                    width:2px;
                    background-color:#FFF;
                    border-right: 1px dashed black;
                    pointer-events: all;                    
                    height:{height+19}px">

        <div style="position:absolute; opacity: 0.5; top:-5px; left:10px; width: 120px">
        <span style="font-size: 13px">{@html math_inline("\\mathcal{R}_t=" + (R0*InterventionAmt).toFixed(2) )}</span> ⟶ 
        </div>

        {#if xScaleTime(InterventionTime) >= 100}
          <div style="position:absolute; opacity: 0.5; top:-2px; left:-97px; width: 120px">
          <span style="font-size: 13px">⟵ {@html math_inline("\\mathcal{R}_0=" + (R0).toFixed(2) )}</span>
          </div>      
        {/if}

      


        </div>
      </div>
     

      <div style="pointer-events: none;
                  position: absolute;
                  top:{height+84}px;
                  left:{0}px;
                  width:{780}px;
                  opacity: 1.0;
                  height:45px;
                  ">
            {#each milestones as milestone}
              <div style="position:absolute; left: {xScaleTime(milestone[0])+8}px; top: -30px;">
                <span style="opacity: 0.3"><Arrow height=30 arrowhead="#circle" dasharray = "2 1"/></span>
                  <div class="tick" style="position: relative; left: 0px; top: 35px; max-width: 130px; color: #BBB; background-color: white; padding-left: 4px; padding-right: 4px">{@html milestone[1]}</div>
              </div>
            {/each}
      </div>       

   </div>

</div>

<div style="height:220px;">
  <div class="minorTitle">
    <div style="margin: 0px 0px 5px 4px" class="minorTitleColumn">Vstupy simulace</div>        
  </div>
  <div class = "row">

    <div class="column">
      <div class="paneltitle">Denní počet COVID pozitivních, kteří dnes zemřeli</div>
      <div class="paneldesc" style="height:60px">V průměru umírá 0.8% lidí, které vir infikuje. Dle počtu úmrtí lze tak zpětně dopočítat počet lidí, kteří virus měli cca 25 dní zpátky. Tento dopočítaný počet pozitivních lidí před 25 dny je pak použit jako počátek simulace, protože tento údaj je spolehlivější, než denní údaje o počtu pozitivních PCR testů</div>
      <div class="slidertext">{startSimulationAtDeaths}</div>
      <input class="range" style="margin-bottom: 8px" type=range bind:value={startSimulationAtDeaths} min=10 max=300 step=1>      
    </div>

    <div class="column">      
      <div class="paneltitle">Reprodukční číslo {@html math_inline("\\mathcal{R}")} do dnešního dne </div>
      <div class="paneldesc" style="height:60px">Toto číslo udává, jak se epidemie vyvíjela během posledních 25 dní. Na základě tohoto údaje a počtu pozitivních před 25 dny (viz vlevo) lze dopočítat, jaká je situace k dnešnímu dni.  </div>      
      <div class="slidertext">{R0}</div>
      <input class="range" type=range bind:value={R0} min=0.01 max=5 step=0.01> 
    </div>    

  </div>
</div>




<div style="position: relative; height: 12px"></div>

<p class = "center">
  <b>Základní reprodukční číslo</b>, někdy také označovano {@html math_inline("\\mathcal{R}_0")} je dáno vlastnostmi viru, a udává kolik dalších lidí v průměru nakazí jeden infekční člověk. Toto základní reprodukční číslo lze však uměle snižovat tím, že viru bráníme v šíření, což přesně modeluje tato kalkulačka. Jak moc které opatření toto reprodukční číslo snižuje je například vidět na skvělé <a href='http://epidemicforecasting.org/containment-calculator'>stránce</a> výzkumného teamu vedeného Janem Kulveitem. 
</p>


<p class = "center">
Tato kalkulačka nepracuje s oficiálními daty ministerstva zdravotnictví (až na údaj o aktuálním denním počtu úmrtí, aby byla schopna určit, v jaké fázi epidemie se nacházíme), ale počítá budoucí vývoj epidemie v ČR pomocí modelu <b><a href="https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model">SEIR</a> </b>(<b>S</b>usceptible → <span style="color:{colors[4]}"><b>E</b></span>xposed → <span style="color:{colors[3]}"><b>I</b></span>nfected → <span><b>R</b></span>emoved) a na základě známého obvyklého průběhu onemocnění v populaci. Tento standardní epidemiologický model pracuje s čtyřmi základními skupinami obyvatel: 

<span>
<ul>
  <li>
    <b>Susceptible</b> (lidé nachylní infekci) jsou lidé, kteří se s virem ještě nepotkali, nemají proti němu imunitu, a mohou se tedy nakazit. Na začátku epidemie jsou v této kategorii všichni lidé v ČR.
  </li>
  <li>
    <b>Exposed</b> (nakažení) jsou lidé, kteří se virem nakazili, ale ještě se u nich vir dostatečně nerozmnožil. Rychlost, s jakou přibývají tito nakažení lidé je dána Reprodukčním číslem {@html math_inline("\\mathcal{R}")}.
  </li>
  <li>
    <b>Infected</b> (infekční) jsou lidé, u kterých se virus již dostatečně rozmnožil a šíří ho dále. Většina lidí pak má i nějaké příznaky, od lekhých po těžší. Malé procento těchto lidí pak vyžaduje nemocniční péči; nicméně i malé procento z velkého počtu může být zásadní.
  </li>
  <li>
    <b>Removed</b> (uzdravení) jsou lidé, kteří nemoc zdárně překonali a mají na nějakou dobu imunitu, čímž brání dalšímu šíření nemoci.
  </li>
</ul>
</span>

Tento model je běžně používán v odborné literatuře [<a href="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30260-9/fulltext">Wu, et. al</a>, <a href = "https://cmmid.github.io/topics/covid19/current-patterns-transmission/wuhan-early-dynamics.html">Kucharski et. al</a>] a jeho dynamika je popsána těmito čtyřmi obyčejnými diferencialními rovnicemi, kde každá popisuje přírustky a úbytky lidí v těchto čtyřech zmíněných skupinách:
<span style="color:#777">{@html ode_eqn}</span>
</p>


<p class = "center">
Jako parametry těchto rovnic jsou použity hodnoty založeny na aktuálním poznání nemoci Covid-19 a jejich typických průběhů. V této kalkulačce byly použity údaje americké <a href='https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html'>CDC</a>, které jsou však v podstatě shodné s údaji udávanými českým <a href='https://www.uzis.cz/res/file/covid/20200917-dusek.pdf'>ÚZIS</a>. Celá kalkulačka je open source, zdrojové kódy jsou volně k dispozici <a href='https://github.com/lukeN86/epcalc'>zde</a>.
</p>



<div style="height:220px;">
  <div class="minorTitle">
    <div style="margin: 0px 0px 5px 4px" class="minorTitleColumn">Dynamika šíření</div>
    <div style="flex: 0 0 20; width:20px"></div>
    <div style="margin: 0px 4px 5px 0px" class="minorTitleColumn">Klinické parametry</div>
  </div>
  <div class = "row">

    <div class="column">
      <div class="paneltitle">Velikost populace</div>
      <div class="paneldesc" style="height:30px"><br></div>
      <div class="slidertext">{format(",")(Math.round(N))}</div>      
    </div>
    <div class="column">
      <div class="paneltitle">Časové parametry</div>
      <div class="paneldesc" style="height:40px">Inkubační doba, {@html math_inline("T_{\\text{inc}}")}.<br></div>
      <div class="slidertext">{(D_incbation).toFixed(2)} dnů</div>      
      <div class="paneldesc" style="height:40px; border-top: 1px solid #EEE; padding-top: 10px">Doba po kterou je člověk infekční, {@html math_inline("T_{\\text{inf}}")}.<br></div>
      <div class="slidertext">{D_infectious} dny</div>      
    </div>
    <div class="column"></div>

    <div style="flex: 0 0 20; width:20px"></div>

    <div class="column">
      <div class="paneltitle">Smrtnost</div>
      <div class="paneldesc" style="height:80px">Infection Fatality Rate, tedy podíl ze všech nakažených lidí (vč. asymptomatických), kteří po nakažení virem zemřou.<br></div>
      <div class="slidertext">{(CFR*100).toFixed(2)} %</div>      
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Doba od projevení prvních příznaků k úmrtí.<br></div>
      <div class="slidertext">{Time_to_death} dnů</div>      
    </div>

    <div class="column">
      <div class="paneltitle">Doba zotavení</div>
      <div class="paneldesc" style="height:30px">Průměrná doba hospitalizace<br></div>
      <div class="slidertext">{D_recovery_severe} dnů</div>      
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Průměrná doba uzdravení doma<br></div>
      <div class="slidertext">{D_recovery_mild} dnů</div>      
    </div>

    <div class="column">
      <div class="paneltitle">Statistiky hospitalizace</div>
      <div class="paneldesc" style="height:30px">Procento nakažených, které vyžaduje hospitalizaci.<br></div>
      <div class="slidertext">{(P_SEVERE*100).toFixed(2)} %</div>      
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Doba od projevení prvních příznaků k hospitalizaci.<br></div>
      <div class="slidertext">{D_hospital_lag} dnů</div>
      
    </div>

  </div>
</div>

<!-- Input data -->
<div style="margin-bottom: 30px">

  <div class="center" style="padding: 10px; margin-top: 3px; width: 925px">
    <div class="legendtext">Link pro sdílení:</div>
    <form>
      <textarea type="textarea" rows="1" cols="5000" style="white-space: nowrap;  overflow: auto; width:100%; text-align: left; height:50px" id="fname" name="fname">{state}</textarea>
    </form>
  </div>
</div>


<p class = "center">
  <b> Autor</b><br>
  Na zakládě kalkulačky, kterou vyvinul <a href='https://gabgoh.github.io/COVID/'>Gabriel Goh</a>, tuto kalkulačku pro situaci v České Republice adaptoval Lukáš Neumann.
</p>

<p class = "center">
  <b> Acknowledgements (English)</b><br>
  Thank you <a href="http://gabgoh.github.io/">Gabriel Goh</a> for the brilliant COVID-19 <a href='https://gabgoh.github.io/COVID/'>calculator</a>! 

  <br/>Original acknowledgements follow:
  <a href = "https://enkimute.github.io/">Steven De Keninck</a> for RK4 Integrator. <a href="https://twitter.com/ch402">Chris Olah</a>, <a href="https://twitter.com/shancarter">Shan Carter
  </a> and <a href="https://twitter.com/ludwigschubert">Ludwig Schubert
  </a> wonderful feedback. <a href="https://twitter.com/NikitaJer">Nikita Jerschov</a> for improving clarity of text. Charie Huang for context and discussion.
  </p>
