// Phase correlation core — same code goes into the page
function makePC(N){
  const LOG = Math.log2(N)|0;
  const rev = new Uint32Array(N);
  for (let i=0;i<N;i++){ let r=0; for(let b=0;b<LOG;b++) r|=((i>>b)&1)<<(LOG-1-b); rev[i]=r; }
  const cosT = new Float32Array(N/2), sinT = new Float32Array(N/2);
  for (let i=0;i<N/2;i++){ cosT[i]=Math.cos(2*Math.PI*i/N); sinT[i]=Math.sin(2*Math.PI*i/N); }
  const hann = new Float32Array(N); for (let i=0;i<N;i++) hann[i]=0.5-0.5*Math.cos(2*Math.PI*i/(N-1));
  const tr = new Float32Array(N), ti = new Float32Array(N);
  function fft1(re, im, off, stride, inv){
    for (let i=0;i<N;i++){ const j=off+rev[i]*stride; tr[i]=re[j]; ti[i]=im[j]; }
    for (let size=2; size<=N; size<<=1){
      const half=size>>1, step=N/size;
      for (let s=0; s<N; s+=size) for (let k=0;k<half;k++){
        const c=cosT[k*step], sn=inv?sinT[k*step]:-sinT[k*step];
        const a=s+k, b=a+half;
        const xr=tr[b]*c-ti[b]*sn, xi=tr[b]*sn+ti[b]*c;
        tr[b]=tr[a]-xr; ti[b]=ti[a]-xi; tr[a]+=xr; ti[a]+=xi;
      }
    }
    for (let i=0;i<N;i++){ const j=off+i*stride; re[j]=tr[i]; im[j]=ti[i]; }
  }
  function fft2(re, im, inv){ for(let y=0;y<N;y++) fft1(re,im,y*N,1,inv); for(let x=0;x<N;x++) fft1(re,im,x,N,inv); }
  // gray: Float32Array N*N (row 0 = top). Fills spectrum (re,im).
  function spectrum(gray, re, im, periodicX){
    let m=0; for (let i=0;i<N*N;i++) m+=gray[i]; m/=N*N;
    for (let y=0;y<N;y++) for (let x=0;x<N;x++){ const i=y*N+x; re[i]=(gray[i]-m)*(periodicX?1:hann[x])*hann[y]; im[i]=0; }
    fft2(re, im, false);
  }
  const cr = new Float32Array(N*N), ci = new Float32Array(N*N);
  // Gaussian low-pass on the whitened cross-spectrum: kills sensor noise that pure phase correlation amplifies
  const lp = new Float32Array(N*N);
  function setLowpass(frac){ const sg = N*frac;
    for (let y=0;y<N;y++){ const fy = y<=N/2?y:y-N; for (let x=0;x<N;x++){ const fx = x<=N/2?x:x-N;
      lp[y*N+x] = Math.exp(-(fx*fx+fy*fy)/(2*sg*sg)); } } }
  setLowpass(0.05);
  // shift of content from A to B (pixels, +x right, +y down) and peak quality 0..1
  function correlate(ar, ai, br, bi){
    for (let i=0;i<N*N;i++){
      const r = ar[i]*br[i] + ai[i]*bi[i];      // conj(A)*B
      const m = ar[i]*bi[i] - ai[i]*br[i];
      const mag = Math.hypot(r, m) + 1e-9;
      cr[i]=r/mag*lp[i]; ci[i]=m/mag*lp[i];
    }
    fft2(cr, ci, true);
    let best=-1e9, bx=0, by=0;
    for (let i=0;i<N*N;i++) if (cr[i]>best){ best=cr[i]; bx=i%N; by=(i/N)|0; }
    const at=(x,y)=>cr[((y+N)%N)*N+((x+N)%N)];
    const sub=(l,c,r)=>{ const d=l-2*c+r; return d<0 ? 0.5*(l-r)/d : 0; };
    let dx = bx + sub(at(bx-1,by), best, at(bx+1,by));
    let dy = by + sub(at(bx,by-1), best, at(bx,by+1));
    if (dx > N/2) dx -= N; if (dy > N/2) dy -= N;
    // peak-to-sidelobe ratio: how far the peak stands above the rest of the surface
    let s1=0, s2=0, n=0; const R=6;
    for (let y=0;y<N;y+=2) for (let x=0;x<N;x+=2){
      let ddx=Math.abs(x-bx); ddx=Math.min(ddx,N-ddx); let ddy=Math.abs(y-by); ddy=Math.min(ddy,N-ddy);
      if (ddx<=R && ddy<=R) continue; const v=cr[y*N+x]; s1+=v; s2+=v*v; n++; }
    const mean=s1/n, sd=Math.sqrt(Math.max(1e-12, s2/n-mean*mean));
    return { dx, dy, peak: best/(N*N), psr: (best-mean)/sd };
  }
  return { spectrum, correlate, setLowpass };
}

// Polar resampling of the amplitude spectrum: rotation of the image becomes a cyclic shift along x.
// Amplitude of a real image is symmetric, so angles cover [0, pi). Rows = radius, cols = angle.
function makePolar(N, P){
  const rmin = 4, rmax = N*0.42;
  const cosA = new Float32Array(P), sinA = new Float32Array(P), rad = new Float32Array(P);
  for (let a=0;a<P;a++){ const t=a*Math.PI/P; cosA[a]=Math.cos(t); sinA[a]=Math.sin(t); }
  for (let j=0;j<P;j++) rad[j] = rmin + j*(rmax-rmin)/(P-1);
  const mag = new Float32Array(N*N);
  function polar(re, im, out){
    for (let i=0;i<N*N;i++) mag[i] = Math.log1p(Math.hypot(re[i], im[i]));
    const at = (x,y) => mag[((y%N+N)%N)*N + ((x%N+N)%N)];
    for (let j=0;j<P;j++){ const r = rad[j], w = r/rmax;          // radius weight = high-pass emphasis
      for (let a=0;a<P;a++){
        const fx = r*cosA[a], fy = r*sinA[a], x0 = Math.floor(fx), y0 = Math.floor(fy), dx = fx-x0, dy = fy-y0;
        const v = at(x0,y0)*(1-dx)*(1-dy) + at(x0+1,y0)*dx*(1-dy) + at(x0,y0+1)*(1-dx)*dy + at(x0+1,y0+1)*dx*dy;
        out[j*P+a] = v*w;
      } }
  }
  return { polar, degPerBin: 180/P };
}
// rotate an N*N gray image by ang (radians, counter-clockwise in a y-up picture) around its centre
function rotateGray(src, dst, N, ang){
  const c = Math.cos(ang), s = Math.sin(ang), h = (N-1)/2;
  for (let y=0;y<N;y++){ const yy = h - y; for (let x=0;x<N;x++){ const xx = x - h;
    // inverse map: dst(p) = src(R^-1 p), y-up coordinates
    const sx = c*xx + s*yy, sy = -s*xx + c*yy;
    const fx = sx + h, fy = h - sy, x0 = Math.floor(fx), y0 = Math.floor(fy);
    if (x0 < 0 || y0 < 0 || x0 >= N-1 || y0 >= N-1){ dst[y*N+x] = NaN; continue; }
    const ax = fx-x0, ay = fy-y0, i = y0*N+x0;
    dst[y*N+x] = src[i]*(1-ax)*(1-ay) + src[i+1]*ax*(1-ay) + src[i+N]*(1-ax)*ay + src[i+N+1]*ax*ay;
  } }
  // fill corners with the mean so they don't add hard edges
  let m=0,n=0; for (let i=0;i<N*N;i++) if (dst[i]===dst[i]){ m+=dst[i]; n++; } m = n? m/n : 0;
  for (let i=0;i<N*N;i++) if (dst[i]!==dst[i]) dst[i]=m;
}
if (typeof module!=='undefined') module.exports = { makePC, makePolar, rotateGray };
