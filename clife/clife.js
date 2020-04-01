
const urlParams = new URLSearchParams(window.location.search);

// CONSTFIG /////////////////////////////////////////////////////////////////////////////////////////

ax = 400;  ay = 300;
az = 1; // number of layers
nT = 2; // number of time instants: 1 = Monte Carlo, 2 = Consistent
mu = 0; // random mutations at birth
breedingType = 1; // 1 = ArithMean, 2 = CrossingOver, 3 = TotalRecombination

gMX = 255; // max integer value for gene (float value = int / gMX = from 0 to 1), storing genes in Uint8 or Uint8Clamped
gMX1 = gMX + 1;
genotypeAccuracy = 100; // genotypes with genes different less than 1/this considered the same

nG = 18; // number of genes

perturn = 1 + round(0.5*ax*ay); // number of point counted per turn for Monte Carlo

// PRESETS /////////////////////////////////////////////////////////////////////////////////////////

if(urlParams.get('preset')) preset = urlParams.get('preset');

     if(preset=='clife')  { ax = 600;  ay = 300;  az = 1;  nT = 2;  mu = 0;  breedingType = 1; }
else if(preset=='xlife')  { ax = 300;  ay = 200;  az = 3;  nT = 2;  mu = 1;  breedingType = 1; }
else if(preset=='mclife') { ax = 800;  ay = 250;  az = 1;  nT = 1;  mu = 0;  breedingType = 1; }

if(urlParams.get('ax')>0) ax = intval(urlParams.get('ax'));
if(urlParams.get('ay')>0) ay = intval(urlParams.get('ay'));
if(urlParams.get('nT')>0) nT = intval(urlParams.get('nT'));
if(urlParams.get('mu')>0) mu = intval(urlParams.get('mu'));

debug = urlParams.get('debug')==1 ? 1 : 0;
stopp = urlParams.get('stopp')==1 ? 1 : 0;
go = urlParams.get('go')==0 ? 0 : 1;

if(debug) { ax = 8;  ay = 4;  stopp = 1;  perturn = 1; }
if(stopp) { go = 0; }

// INIT RULES /////////////////////////////////////////////////////////////////////////////////////////

function initRules() {
  rules = genlims = coarses = [];
  if(urlParams.get('random')==1) { // random rules to look for something new and interesting
    var rule0 = [];
    for(var v=0; v<7; v++) {
      rule0[v] = randbs(true) + ':' + randbs(false);
    }
    rules[0] = rule0;
    console.log(rule0);
  }
  else if(nT==2 && az==1) { // classic Life-like (clife)
    rules = [
      [
        '37:23',  // DryLife
        //'3:023',  // DotLife
        '36:125',  // 2x2
        '3:23',  // Conway's Life
        '36:23',  // HighLife
        '357:238',  // Pseudo Life
        '38:23',  // Pedestrian Life
        '368:238',  // LowDeath
        '38:238',  // HoneyLife
        '3:238',  // EightLife
        
        //'357:1358',  // Amoeba
        //'35678:5678',  // Diamoeba
        //'34:456',  // Bacteria
        //'3:45678',  // Coral
        //'34578:456',  // Gems Minor
        //'36:235',  // Blinker Life
      ]
    ];
  }
  else if(nT==2) { // Multi-layer (xlife)
    rules = [
      [ // plants
        //  '357:1358',  // Amoeba
        '35678:5678',  // Diamoeba
        //, '34:456',  // Bacteria
        //, '3:45678',  // Coral
        //, '34578:456',  // Gems Minor
        //, '36:235',  // Blinker Life
      ],
      [ // herbivores
        '37:235',
        '368:237',
        '358:235',
      ],
      [ // carnivores
        '2:013',
        '12:02',
      ]
    ];
    //rules[1] = ['245:235', '', '28:23', '', '248:37', '']; // fire!
    //rules[1] = ['25:7', '']; // diam acid
    rules[1] = ["247:06", "6:015", "2:2", "36:3"];
    if(urlParams.get('random')==2) {
      var rule0 = [];
      for(var v=0; v<7; v++) {
        rule0[v] = randbs(true) + ':' + randbs(false);
      }
      rules[1] = rule0;
      console.log(rule0);
    }
    genlims = [10*gMX, 7*gMX, 5*gMX]; // limit on sum of genes in genotype
    coarses = [2, 1, 2]; // rounds genes: spectrum of possible gene values = 256 / 2^coarse
  }
  else if(nT==1 && az==1) { // Monte Carlo (mclife)
    rules = [
      [
        '348:247',
        '35:348',
        '3:1348',
        '34:24',
        
        '467:4568',
        '5:456',
        '568:345',
        
        '26:38',
        '27:07',
        '27:38',
        '2:3',
        '258:04',
        //'28:057',
      ]
    ];
  }
}

function initialFill() {
  var lx = 1.0, ly = 0.8;
  
  for(var z=0; z<az; z++) {
    for(var x=round(ax/2-ax*lx/2); x<round(ax/2+ax*lx/2); x++) {
      for(var y=round(ay/2-ay*ly/2); y<round(ay/2+ay*ly/2); y++) {
        if(z==2 && y<ay/2) continue; // not too much predators
        var density = round((1 - Math.abs(2*y/ay-1)/ly)*10)/10;  // Math.sin(Math.PI * x / ax);
        if(Math.random()<=density) {
          var r = floor(rules[z].length * (x-round(ax/2-ax*lx/2)) / ax / lx);  // r = rnd0(rules[z].length);
          if(!rules[z][r]) continue;
          BornCell(z, x, y, bs2arr(rules[z][r]), true);
        }
      }
    }
  }
}

// GAME MECHANICS /////////////////////////////////////////////////////////////////////////////////////////

function BornCell(z, x, y, genes, inc=true) {
  genes = new Uint8ClampedArray(genes);
  var hsh = arr2hsh(genes, coarses[z]);
  if(!G[hsh]) {
    G[hsh] = genes;
    Gcount ++;
  }
  
  F[C(T1,z,x,y)] = hsh;
  
  if(inc) IncNB(z, x, y);
}

function DieCell(z, x, y, inc=true) {
  F[C(T1,z,x,y)] = 0;
  
  if(inc) IncNB(z, x, y, -1);
}

function IncNB(z, x, y, d=1) {
  var dx, dy, x0, y0;
  for(dx=-1; dx<=1; dx++) {
    x0 = x + dx;  if(x0<0) x0 = ax - 1;  if(x0>=ax) x0 = 0;
    for(dy=-1; dy<=1; dy++) {
      if(dx==0 && dy==0) continue; // you yourself are not your neighbor
      y0 = y + dy;  if(y0<0) y0 = ay - 1;  if(y0>=ay) y0 = 0;
      NB[C(T1,z,x0,y0)] += d;
    }
  }
}

function MutateGenes(genes, z=0) {
  var r = Math.random();
  var num = 0;
       if(r<0.001) num = 3;
  else if(r<0.01)  num = 2;
  else if(r<0.1)   num = 1;
  else return genes;
  
  if(genlims[z]) {
    var sum = 0;
    for(var i=0; i<nG; i++) sum += genes[i];
  }
  
  for(var n=0; n<num; n++) {
    var gn = rnd0(nG); // choosing random gene
    if(gn<=1) return genes; // genes 0 and 1 should never be positive
    var newgene = MutateGene(genes[gn]);
    if(genlims[z]) {
      if(sum + newgene - genes[gn] > genlims[z]) return genes; // reached limit on gene sum in genotype
      sum += newgene - genes[gn];
    }
    genes[gn] = newgene;
  }
  
  return genes;
}

function MutateGene(g) {
  var r = Math.random();
  var sgn;
  if(r<=0.1) { sgn = -1; }
  else if(r>=0.9) { sgn = 1;   r = 1 - r; }
  else return g;
  g = g - sgn * Math.log10(r) * 0.1 * gMX; // 10% of +0.1, 1% of +0.2, etc
  if(g<0) g = 0;
  else if(g>gMX) g = gMX;
  return gMX>100 ? round(g) : g;
}

function calcCell(z, x, y) {
  var dx, dy, i, j, r;
  
  var nneib = NB[C(T0,z,x,y)];
  
  var FCT0zxy = F[C(T0,z,x,y)];
  if(FCT0zxy) { // living cell - will it survive?
    var psurv = G[FCT0zxy][9+nneib]; // probability to survive
    
    var surv; // will it survive?
    if(!psurv) surv = false;
    else if(psurv>=gMX) surv = true;
    else surv = rnd0(gMX) < psurv ? true : false;
    
    if(surv) {
      if(nT>1) {
        F[C(T1,z,x,y)] = FCT0zxy;
        IncNB(z, x, y);
      }
      if(z>0 && F[C(T1,z-1,x,y)]) { // animal eats some grass when survives
        r = Math.random();
        if(z==1 && r<0.005 || z>=2 && r<0.05) {
          DieCell(z-1, x, y, true);
        }
      }
    }
    else {
      DieCell(z, x, y,
        nT==1 ? true : false // if nT==2 we start NB from fill(0)
      );
    }
  }
  else if(!nneib) { // dead and no neighbors -> guaranteed dead
    if(nT>1) DieCell(z, x, y, false);
  }
  else if(nneib<2 && z<2) { // dead and neighbors<2 and not a predator -> stay dead
    if(nT>1) DieCell(z, x, y, false);
  }
  else if(z>0 && !F[C(T1,z-1,x,y)]) { // in X-life for cell to born we need previous layer cell to be alive (next layer cell will eat it)
    if(nT>1) DieCell(z, x, y, false);
  }
  else { // dead cell - will it be born?
    // for a new cell rules are inherited from neighbors
    // finding gene responsible for birth when neighbors count = nneib
    var neibgenes = [], nn = 0;
    var pborn = 0; // probability to born
    for(dx=-1; dx<=1; dx++) {
      x0 = x + dx;  if(x0<0) x0 = ax - 1;  if(x0>=ax) x0 = 0;
      for(dy=-1; dy<=1; dy++) {
        if(dx==0 && dy==0) continue;
        y0 = y + dy;  if(y0<0) y0 = ay - 1;  if(y0>=ay) y0 = 0;
        var f = F[C(T0,z,x0,y0)];
        if(f) {
          neibgenes[nn] = G[f];
          pborn += neibgenes[nn][nneib];
          nn ++;
        }
      }
    }
    pborn = pborn / nneib;
    
    var born; // will it be born?
    if(!pborn) born = false;
    else if(pborn>=gMX) born = true;
    else born = rnd0(gMX) < pborn ? true : false;
    
    if(born) {
      // breeding - defining child's genes from parents' genes
      var genes = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
      if(breedingType==1) { // ArithMean
        for(i=0; i<nneib; i++) {
          for(j=0; j<nG; j++) {
            genes[j] += neibgenes[i][j] / nneib;
          }
        }
      }
      else if(breedingType==2) { // CrossingOver
        var cg = 0; // current child's gene number, will increase from 0 up to nG
        for(i=0; i<nG; i++) {
          var tp = rnd0(nneib); // parent
          var tl = rnd(1, nG); // length of this parent's gene sequence inherited by child
          for(j=0; j<tl; j++) {
            genes[cg] = neibgenes[tp][cg];
            cg++;  if(cg>=nG) break;
          }
          if(cg>=nG) break;
        }
      }
      else if(breedingType==3) { // TotalRecombination
        // each gene selects from a random parent
      }
      
      if(mu) genes = MutateGenes(genes, z);
      
      BornCell(z, x, y, genes, true);
      
      if(z==1 || (z>=2 && Math.random()<0.5)) { // animal eats grass to be born (condition for grass existence is above)
        DieCell(z-1, x, y, true);
      }
    }
    else {
      if(nT>1) DieCell(z, x, y, false);
    }
  }
  
}

function calcWorld() {
  if(!nturn) ShowStats();
  
  drawWorld();
  
  if(!go) return 0;
  
  if(stopp) go = 0;
  
  if(nT>1) {
    if(T1==1) { T1 = 0;  T0 = 1; } else { T1 = 1;  T0 = 0; } // switching between previous and current moment fields
    NB.fill(0, Mt*T1, Mt*T1+Mt);
  }
  
  var z, x, y;
  if(nT==1) {
    for(var l=0; l<perturn; l++) {
      z = rnd0(az);
      x = rnd0(ax);
      y = rnd0(ay);
      
      calcCell(z, x, y);
    }
  }
  else {
    for(z=0; z<az; z++) {
      for(x=0; x<ax; x++) {
        for(y=0; y<ay; y++) {
          calcCell(z, x, y);
        }
      }
    }
  }
  
  nturn++;
  
  // G cleanup every X turns
  if(nturn>0 && !(nturn%10)) {
    var GN = {};
    for(z=0; z<az; z++) {
      for(x=0; x<ax; x++) {
        for(y=0; y<ay; y++) {
          var hsh = F[C(T1,z,x,y)];
          if(!hsh) continue;
          if(!GN[hsh]) GN[hsh] = 0;  GN[hsh] ++; // counting cells with that genotype to delete unused from G
        }
      }
    }
    
    for(var hsh in G) {
      if(GN[hsh]>0) continue;
      delete G[hsh];
      Gcount --;
    }
  }
  
  // calc fps
  var date1 = new Date;
  var timer = date1 - date0;  date0 = date1;
  
  timerv += timer;  timern++;
  if(timern>=50 || stopp) {
    ShowStats();
  }
  
  var delay = 15 - timer;  if(delay<0) delay = 0;
  if(go) setTimeout(calcWorld, delay);
}

// STATISTICS /////////////////////////////////////////////////////////////////////////////////////////

tracked = [];  prevpoints = [];  infostep = -1;
function ShowStats() {
  var s1 = s2 = '';
  var z, x, y, i, j, d, idx, zidx;
  
  infostep ++;
  
  // show fps
  var ms = timerv / timern;
  timerv = 0;  timern = 0;
  s1 += 'turn = ' + nturn + '<br>';
  s1 += 'fps = ' + round(1000/ms) + '<br>';
  
  // calc species stats
  var specstat = [], ttl = 0;  for(z=0; z<az; z++) specstat[z] = [];
  for(z=0; z<az; z++) {
    for(x=0; x<ax; x++) {
      for(y=0; y<ay; y++) {
        if(!F[C(T1,z,x,y)]) continue;
        var ga = G[F[C(T1,z,x,y)]];
        if(!ga) continue; // genotype already deleted from list
        idx = '';
        for(i=0; i<nG; i++) {
          var gi = ga[i] / gMX;
               if(gi<0.125) d = '0';
          else if(gi<0.375) d = 'a';
          else if(gi<0.625) d = 'b';
          else if(gi<0.875) d = 'c';
          else              d = '1';
          idx += (idx?' ':'') + d;
        }
        if(!specstat[z][idx]) specstat[z][idx] = 0;
        specstat[z][idx] ++;
        ttl ++;
      }
    }
  }
  s1 += 'live cells = ' + ttl + '<br>';
  
  if(0) {
    console.log(Object.keys(G).length);
    //console.log(G);
    console.log(GN);
  }
  
  s1 += 'genotypes = ' + Gcount + '<br>';
  
  // tracking all species ever reached top-3 in their z-plane or filled 0.1% of field's area
  for(z=0; z<az; z++) {
    var sorted = arsort_keys(specstat[z]);
    var l = sorted.length, lmt = round(0.001*ax*ay);
    for(i=0; i<l; i++) {
      idx = sorted[i];  if(!idx) break;  if(!specstat[z][idx]) break;
      if(i<3 || specstat[z][idx]>lmt) {
        zidx = z + ': ' + idx;
        tracked[zidx] = specstat[z][idx];
      }
      else break;
    }
  }
  
  // updating all tracked values
  for(zidx in tracked) {
    z = zidx.substring(0, 1);  idx = zidx.substring(3);
    if(specstat[z][idx]) tracked[zidx] = specstat[z][idx];
    else tracked[zidx] = 0;
  }
  
  var trsorted = arsort_keys(tracked);
  
  s2 += 'dominating genotypes:<br>';
  for(i in trsorted) {
    zidx = trsorted[i];
    z = zidx.substring(0, 1);  idx = zidx.substring(3);
    var genes = [], g = 0;
    var chars = z + ':', c = '';
    var splt = idx.split(' ');
    for(j=0; j<nG; j++) {
           if(splt[j]=='0') { g = 0;     c = '0'; }
      else if(splt[j]=='a') { g = 0.25;  c = '&frac14;'; }
      else if(splt[j]=='b') { g = 0.5;   c = '&frac12;'; }
      else if(splt[j]=='c') { g = 0.75;  c = '&frac34;'; }
      else if(splt[j]=='1') { g = 1;     c = '1'; }
      else                  { g = 0;     c = '?'; }
      genes[j] = g * gMX;  chars += ' ' + c ;
    }
    var clr = ColorifyGenes(genes, z);
    s2 +=
      '<span style="color:rgb('+clr.r+','+clr.g+','+clr.b+'); background:#000;">' + chars + '</span>'
      + ' = ' + tracked[zidx]
      + ' = ' + round(tracked[zidx]/ttl*100) + '%'
      + ' &nbsp; ' + arr2bs(genes)
      + '<br>';
    
    if(infostep<zoom*ax) {
      var xx = infostep;
      var yy = scnv_height - round(Math.log2(tracked[zidx]) / Math.log2(ax*ay) * scnv_height); // Math.log2 or cbrt here
      var style = 'rgb('+clr.r+','+clr.g+','+clr.b+')';
      if(prevpoints[zidx] && xx>0) {
        sctxs[z].beginPath();
        sctxs[z].strokeStyle = style;
        sctxs[z].moveTo(xx-1, prevpoints[zidx]);
        sctxs[z].lineTo(xx, yy);
        sctxs[z].stroke();
      }
      else {
        sctxs[z].fillStyle = style;
        sctxs[z].fillRect(xx, yy, 1, 1);
      }
      prevpoints[zidx] = yy;
    }
    
    if(!tracked[zidx]) {
      delete tracked[zidx];
    }
  }
  
  document.getElementById('stxt1').innerHTML = s1;
  document.getElementById('stxt2').innerHTML = s2;
}

// RULES /////////////////////////////////////////////////////////////////////////////////////////

function randbs(no012=false) {
  var r = '', d0 = -1;
  for(var i=0; i<9; i++) {
    var d = rnd(no012?2:0, 9);
    if(d<=d0) break;
    d0 = d;
    r += d.toString();
  }
  return r;
}

function bs2arr(bs) {
  var nb, ns, ab, as, born, surv, k;
  [nb, ns] = bs.split(':');
  ab = nb.split('');  born = [0, 0, 0, 0, 0, 0, 0, 0, 0];  if(ab) for(k in ab) born[ab[k]] = gMX;
  as = ns.split('');  surv = [0, 0, 0, 0, 0, 0, 0, 0, 0];  if(as) for(k in as) surv[as[k]] = gMX;
  return [...born, ...surv];
}

function arr2bs(arr) {
  var bs = '';
  for(var i=0; i<nG; i++) {
    if(i==9) bs += ':';
    if(arr[i]>=0.99*gMX) bs += i%9;
    else if(arr[i]<=0.01*gMX) bs += '';
    else return '';
  }
  return bs;
}

// MATHS /////////////////////////////////////////////////////////////////////////////////////////

function str2hsh(s) {
  var hash = 0, i, l, ch;
  for(i=0,l=s.length; i<l; i++) {
    ch = s.charCodeAt(i);
    hash = ((hash << 5) - hash) + ch;
    hash |= 0;
  }
  return hash;
}
function arr2hsh(ar, coarse=2) {
  var hash = 0, i, l, ch;
  for(i=0,l=ar.length; i<l; i++) {
    //ch = round(ar[i]*genotypeAccuracy);  //if(!ch.isInteger) ch = float2IEEE(ch);
    ch = ar[i] >> coarse;
    hash = ((hash << 5) - hash) + ch;
    hash |= 0;
  }
  return hash;
}
function float2IEEE(x) {
  var buf = new ArrayBuffer(4);
  (new Float32Array(buf))[0] = x;
  return (new Uint32Array(buf))[0];
}

function sqr2(x) { return x*x; }
function floor(x){ return Math.floor(x); }
function round(x){ return Math.round(x); }
function round4(x) { return Math.round(x * 10000 + Number.EPSILON) / 10000; }
function intval(x) { return parseInt(x, 10); }
function rnd0(n) { return Math.floor(n*Math.random()); }
function rnd(a, b) { return Math.floor(Math.random()*(b-a)) + a; }

function arsort_keys(obj) {
  var keys = Object.keys(obj);
  return keys.sort(function(a,b){return obj[b]-obj[a]});
}

sqr255 = 255 * 255;

// GRAPHICS /////////////////////////////////////////////////////////////////////////////////////////

function ColorifyGenes(genes, z) {
  var r, g, b, l;
  r = (genes[4] + genes[8]) / 2 + (genes[ 9] + genes[12] + genes[17]) / 3;  // if(z==2) r += 0.5;
  g = (genes[5] + genes[7]) / 2 + (genes[10] + genes[14] + genes[16]) / 3;  // if(z==0) g += 0.5;
  b = (genes[6] + genes[6]) / 2 + (genes[11] + genes[13] + genes[15]) / 3;  // if(z==1) b += 0.5;
  l = Math.sqrt(r*r + g*g + b*b);
  if(l) {
    r = 255 * r / l;
    g = 255 * g / l;
    b = 255 * b / l;
  }
  else { r = g = b = 0; }
  return {"r": r, "g": g, "b": b};  // no round needed for Uint8ClampedArray?
}

function GetPixel(x, y) {
  var s = 4 * (y * axzoom + x);
  var r = fimg.data[s+0];
  var g = fimg.data[s+1];
  var b = fimg.data[s+2];
  return [r, g, b];
}

function SetPixel(x, y, r, g, b) {
  if(isLittleEndian) {
    fbufU[y * axzoom + x] =
      (255 << 24) |
      (b   << 16) |
      (g   <<  8) |
       r;
  }
  else {
    fbufU[y * axzoom + x] =
      (r   << 24) |
      (g   << 16) |
      (b   <<  8) |
       255;
  }
  /*
    // SetPixel without ArrayBuffer if it is not supported
    var s = 4 * (y * axzoom + x);
    fimg.data[s+0] = r;
    fimg.data[s+1] = g;
    fimg.data[s+2] = b;
    fimg.data[s+3] = 255;
  */
}

function DrawCell(z, x, y, r, g, b, f=false) {
  var mode, xx, yy;
  
  if(zoom==1) { // cell = 1 pixel
    mode = 1;
    xx = x;  yy = y;
  }
  else if(zoom==2 && az==1) { // cell = 2*2 square
    mode = 2;
    xx = zoom * x;  yy = zoom * y;
  }
  else if(zoom==2) { // cell = mixture of pixels from all z-layers
    mode = 3;
    xx = zoom * x;  yy = zoom * y;
         if(z==0) {  }
    else if(z==1) { xx += 1; }
    else if(z==2) { xx += 1;  yy += 1; }
    else if(z==3) { yy += 1; }
  }
  else { // cell = big square (mainly for debug version)
    mode = 4;
    xx = zoom * x;  yy = zoom * y;
  }
  
  if(f) { // fading pixels for dying cells
    [r, g, b] = GetPixel(xx, yy);
    var rgb = r + g + b; // better to use sqrt(r*r+g*g+b*b), but its slower
    if(rgb<60) { // totally decayed pixel
      r = g = b = 0;
    }
    else if(rgb>250) { // have just died
      r = (0.5*r);
      g = (0.5*g);
      b = (0.5*b);
    }
    else { // decaying pixel
      r = (0.9*r);
      g = (0.9*g);
      b = (0.9*b);
    }
  }
  
  if(mode==1) {
    SetPixel(x, y, r, g, b);
  }
  else if(mode==2) {
    SetPixel(xx,   yy,   r, g, b);
    SetPixel(xx+1, yy,   r, g, b);
    SetPixel(xx,   yy+1, r, g, b);
    SetPixel(xx+1, yy+1, r, g, b);
  }
  else if(mode==3) {
    SetPixel(xx, yy, r, g, b);
         if(z==0) {  }
    else if(z==1) { SetPixel(xx-1, yy+1, r, g, b); }
    else if(z==2) { if(r>0||g>0||b>0) SetPixel(xx-1, yy, r, g, b); }
    else if(z==3) {  }
  }
  else {
    for(var i=0; i<zoom; i++) {
      xx = zoom * x + i;
      for(var j=0; j<zoom; j++) {
        yy = zoom * y + j;
        SetPixel(xx, yy, r, g, b);
      }
    }
  }
}

function drawWorld() {
  var x, y, color, r, g, b;
  for(z=0; z<az; z++) {
    for(x=0; x<ax; x++) {
      for(y=0; y<ay; y++) {
        var fc = F[C(T1,z,x,y)];
        if(fc) {
          var gf = G[fc];
          color = ColorifyGenes(gf, z);
          DrawCell(z, x, y, color.r, color.g, color.b);
        }
        else {
          DrawCell(z, x, y, 0, 0, 0, true);
        }
        
      }
    }
  }
  fimg.data.set(fbuf8);
  fctx.putImageData(fimg, 0, 0);
}

// INIT WORLD /////////////////////////////////////////////////////////////////////////////////////////

function initVar() {
  nturn = 0;
  date0 = new Date;
  timerv = timern = 0;
  
  T0 = 0;
  T1 = nT>1 ? 1 : 0; // nT==2 means: we keep two copies of The Field, T - current moment, T0 - previous moment
  
  F  = new Int32Array(nT * az * ax * ay); // The Field
  NB = new Uint8Array(nT * az * ax * ay); // neighbors count
  
  Mt = az * ax * ay;
  Mz = ax * ay;
  Mx = ay;
  My = 1; // multipliers to map coords [t][z][x][y] to 1D array index
  C = (t, z, x, y) => t*Mt + z*Mz + x*Mx + y;
  
  G = {}; // genotypes list
  Gcount = 0; // counters for G
}



// INTERFACE INIT /////////////////////////////////////////////////////////////////////////////////////////

function initClicks() {
  fcnv.addEventListener('click', function(ev) {
    go = false;
    var x = floor(ev.offsetX / zoom);
    var y = floor(ev.offsetY / zoom);
    for(var z=0; z<az; z++) {
      console.log(F[C(T1,z,x,y)]);
    }
    SetPixel(x*zoom, y*zoom, 255, 255, 255);  fctx.putImageData(fimg, 0, 0);
  });
}

function initDraw() {
  if(document.body.clientWidth < document.body.clientHeight) [ax, ay] = [ay, ax];
  var zoomx = floor(0.9 * document.body.clientWidth  / ax);  if(zoomx<1) zoomx = 1;
  var zoomy = floor(0.9 * document.body.clientHeight / ay);  if(zoomy<1) zoomy = 1;
  zoom = Math.min(zoomx, zoomy);
  if(az>1) zoom = 2;
  
  // field canvas
  fcnv = document.getElementById('fieldcnv');
  fcnv.width  = zoom * ax;  fcnv.style.width  = fcnv.width  + 'px';
  fcnv.height = zoom * ay;  fcnv.style.height = fcnv.height + 'px';
  fctx = fcnv.getContext('2d');
  fimg = fctx.createImageData(fcnv.width, fcnv.height);
  // array buffer for faster SetPixel (one 32-bit is better than 4 operations on 8-bit)
  // https://hacks.mozilla.org/2011/12/faster-canvas-pixel-manipulation-with-typed-arrays/
  fbuf  = new ArrayBuffer(fimg.data.length); 
  fbuf8 = new Uint8ClampedArray(fbuf);
  fbufU = new Uint32Array(fbuf);
  // determine whether Uint32 is little- or big-endian
  isLittleEndian = true;
  fbufU[1] = 0x0a0b0c0d;
  if(fbuf[4]===0x0a && fbuf[5]===0x0b && fbuf[6]===0x0c && fbuf[7]===0x0d) isLittleEndian = false;
  
  axzoom = ax * zoom; // often used in SetPixel
  
  // stats canvas
  scnv_height = 200;
  scnvs = sctxs = [];
  sdiv = document.getElementById('statcanvas');
  sdiv.innerHTML = 'species population (log scale):<br>';
  for(var z=0; z<az; z++) {
    scnvs[z] = document.createElement('canvas');
    sdiv.appendChild(scnvs[z]);
    scnvs[z].width  = zoom * ax;    scnvs[z].style.width  = scnvs[z].width  + 'px';
    scnvs[z].height = scnv_height;  scnvs[z].style.height = scnvs[z].height + 'px';
    scnvs[z].style.margin = '0 0 5px 0';
    sctxs[z] = scnvs[z].getContext('2d');
  }
  
  document.getElementById('pausecont').style.width = zoom * ax + 'px';
}

// START IT! /////////////////////////////////////////////////////////////////////////////////////////

function init() {
  initDraw();
  initClicks();
  initVar();
  initRules();
  initialFill();
  setTimeout(calcWorld, 10);
}

window.onload = init;

///////////////////////////////////////////////////////////////////////////////////////////
