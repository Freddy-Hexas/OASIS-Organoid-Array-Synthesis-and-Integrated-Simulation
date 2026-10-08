/* Result-only assessment. These controls never enter a routing request. */
(() => {
  'use strict';
  const byId = id => document.getElementById(id);
  const ids = ['assessmentRadiusMm', 'assessmentCenterXmm', 'assessmentCenterYmm', 'assessmentDistanceUm', 'assessmentResolution'];
  const settingsStorage = 'gds-fixed-assessment-v1', comparisonStorage = 'gds-assessment-comparison-v1';
  let job = null, detail = null, report = null, lastKey = null, epoch = 0, controller = null, dirty = false;
  let comparisons = [], field = null, canvasBounds = null;
  const number = (x, digits = 2) => x == null ? '—' : Number(x).toLocaleString('zh-CN', { maximumFractionDigits: digits });
  const percent = x => `${number(100*x, 2)}%`;
  const status = (message, error = false) => { byId('distributionStatus').textContent = message; byId('distributionStatus').classList.toggle('error', error); };
  const node = (tag, text, className) => { const n = document.createElement(tag); if (text != null) n.textContent = text; if (className) n.className = className; return n; };
  function readSettings() {
    const values = ids.map(id => { const input = byId(id); if (!input.value.trim()) throw new Error('请填写全部评价参数'); return Number(input.value); });
    if (!values.every(Number.isFinite) || values[0] <= 0 || values[3] < 0 || ![192,384,768].includes(values[4])) throw new Error('目标半径须大于 0；覆盖距离可为 0；圆心须为有限数值');
    const target = { radius_um: values[0]*1000, center_um: [values[1]*1000, values[2]*1000], coverage_distance_um: values[3], resolution: values[4] };
    if (![target.radius_um, ...target.center_um].every(Number.isFinite)) throw new Error('评价尺寸超出数值范围');
    return target;
  }
  const targetKey = t => JSON.stringify([t.radius_um, ...t.center_um, t.coverage_distance_um, t.resolution]);
  function saveSettings() { try { localStorage.setItem(settingsStorage, JSON.stringify(ids.map(id => byId(id).value))); } catch {} }
  try {
    const previous = JSON.parse(localStorage.getItem(settingsStorage));
    if (Array.isArray(previous) && previous.length === ids.length) {
      const initial = ids.map(id => byId(id).value);
      previous.forEach((value, i) => { byId(ids[i]).value = String(value); });
      try { readSettings(); } catch { initial.forEach((value, i) => { byId(ids[i]).value = value; }); }
    }
    const snapshots = JSON.parse(localStorage.getItem(comparisonStorage));
    if (Array.isArray(snapshots)) comparisons = snapshots.filter(s => s && s.target && s.job_id && s.nearest_neighbor && s.maximum_uncovered && s.coverage && s.curve).slice(-8);
  } catch {}

  function invalidate() {
    dirty = true; lastKey = null; ++epoch; controller?.abort(); report = null;
    byId('distributionContent').hidden = true;
    byId('distributionBadge').textContent = '参数待应用';
    status('评价参数已改变。点击“计算 / 更新评价”后使用新的固定目标区域；已有布线不变。');
  }
  ids.forEach(id => byId(id).addEventListener('input', invalidate));
  byId('assessmentUseCenter').addEventListener('click', () => {
    const center = detail?.geometry?.electrode_region_reference_um;
    if (!Array.isArray(center)) return;
    byId('assessmentCenterXmm').value = String(center[0]/1000);
    byId('assessmentCenterYmm').value = String(center[1]/1000);
    invalidate();
  });
  async function calculate(force = false) {
    if (job?.status !== 'complete' || !job.result) return;
    let target;
    try { target = readSettings(); } catch (error) { status(error.message, true); return; }
    const key = `${job.id}:${targetKey(target)}`;
    if (!force && key === lastKey) return;
    lastKey = key; const token = ++epoch; controller?.abort(); controller = new AbortController();
    report = null; byId('distributionContent').hidden = true; dirty = false; saveSettings();
    byId('distributionBadge').textContent = '计算中…';
    status('正在评价最终电极坐标；只读取已有结果，不启动布线任务。');
    const params = new URLSearchParams({ radius_um: target.radius_um, center_x_um: target.center_um[0], center_y_um: target.center_um[1], coverage_distance_um: target.coverage_distance_um, resolution: target.resolution });
    try {
      const response = await fetch(`/api/distribution/${job.id}?${params}`, {cache:'no-store', signal: AbortSignal.any([controller.signal, AbortSignal.timeout(90000)])});
      const data = await response.json();
      if (!response.ok) throw new Error(data.error || `HTTP ${response.status}`);
      if (token !== epoch) return;
      report = data; renderReport();
      byId('distributionBadge').textContent = '已有布局 · 评价完成';
      status('结果已更新。所有指标使用下方显示的同一个目标圆；几何覆盖距离不等同于已校准的生物信号探测距离。');
    } catch (error) {
      if (token !== epoch) return;
      byId('distributionBadge').textContent = '评价未完成';
      status(`分布评价失败：${error.message}。可点击重新计算。`, true);
    }
  }
  byId('assessmentApply').addEventListener('click', () => calculate(true));
  function renderReport() {
    const t = report.target, nn = report.nearest_neighbor, h = report.maximum_uncovered, coverage = report.coverage;
    byId('distributionContent').hidden = false;
    byId('distributionTarget').textContent = `固定目标圆：圆心 (${number(t.center_um[0]/1000,6)}, ${number(t.center_um[1]/1000,6)}) mm · 半径 ${number(t.radius_um/1000,6)} mm · 覆盖距离 ${number(t.coverage_distance_um,6)} μm · ${t.resolution} × ${t.resolution}。`;
    byId('distributionCount').textContent = number(report.electrode_count,0);
    byId('distributionInside').textContent = `圆内 ${report.inside_target_count} 个 · 圆外 ${report.electrode_count-report.inside_target_count} 个`;
    byId('distributionCV').textContent = nn.cv == null ? '未定义' : number(nn.cv,4);
    byId('distributionNNDetail').textContent = nn.cv == null ? '至少需要 2 个电极且平均间距大于 0' : `全部电极：均值 ${number(nn.mean_um)} μm · 标准差 ${number(nn.std_um)} μm`;
    byId('distributionMaximum').textContent = h.unbounded ? '∞（无电极）' : `≈ ${number(h.estimate_um)} μm`;
    byId('distributionMaximumDetail').textContent = h.unbounded ? '目标区域没有任何可供计算距离的电极' : `区间 ${number(h.lower_um)}–${number(h.upper_um)} μm`;
    byId('distributionCoverageLabel').textContent = `${number(t.coverage_distance_um)} μm 距离覆盖率`;
    byId('distributionCoverage').textContent = `≈ ${percent(coverage.estimate)}`;
    byId('distributionCoverageDetail').textContent = `区间 ${percent(coverage.lower)}–${percent(coverage.upper)}`;
    byId('distributionAccuracy').textContent = report.heatmap ?
      `网格步长 ${number(report.numerics.cell_step_um,3)} μm，距离空间误差不超过约 ${number(report.numerics.distance_uncertainty_um,3)} μm；圆与网格交集用解析面积加权。区间包含空间采样误差，仍使用浮点计算。坐标来自已通过导出回读的最终报告；统计全部 ${report.electrode_count} 个电极，圆外点也可贡献圆内覆盖。此处评价二维欧氏距离，三维成形及组织接触需另行评价。` :
      '没有电极：最近邻指标未定义，最大未覆盖距离为无穷，所有有限距离的覆盖率为 0。';
    renderComparisons(); drawHeatmap();
  }
  const matchingComparisons = () => report ? comparisons.filter(item => targetKey(item.target) === targetKey(report.target)) : [];
  function renderComparisons() {
    const matching = matchingComparisons(), rows = byId('distributionComparisonRows'); rows.replaceChildren();
    byId('distributionComparisonNote').textContent = `比较固定圆心、半径、覆盖距离及计算精度，并同时报告电极数量。当前匹配 ${matching.length} 个保存项${comparisons.length > matching.length ? `；${comparisons.length-matching.length} 个区域或尺度不同的保存项已排除` : ''}。`;
    if (!matching.length) { const tr=node('tr'),td=node('td','点击“加入本组比较”，保存当前布局的评价快照。');td.colSpan=6;tr.append(td);rows.append(tr); }
    for (const item of matching) {
      const tr=node('tr'), title=node('td',item.input_name || item.job_id); title.append(node('small',`任务 ${item.job_id}`));tr.append(title);
      [String(item.electrode_count),String(item.inside_target_count),number(item.nearest_neighbor.cv,4),item.maximum_uncovered.unbounded ? '∞' : `${number(item.maximum_uncovered.lower_um)}–${number(item.maximum_uncovered.upper_um)} μm`,`${percent(item.coverage.estimate)} [${percent(item.coverage.lower)}–${percent(item.coverage.upper)}]`].forEach(text=>tr.append(node('td',text)));
      rows.append(tr);
    }
    drawCurve();
  }
  byId('assessmentAddComparison').addEventListener('click', () => {
    if (!report) return;
    const snapshot = Object.fromEntries(['job_id','input_name','input_sha256','target','electrode_count','inside_target_count','nearest_neighbor','maximum_uncovered','coverage','curve'].map(key=>[key,report[key]]));
    comparisons = comparisons.filter(s=>s.job_id!==snapshot.job_id || targetKey(s.target)!==targetKey(snapshot.target));
    comparisons.push(snapshot); comparisons=comparisons.slice(-8);
    try { localStorage.setItem(comparisonStorage,JSON.stringify(comparisons)); } catch {}
    renderComparisons();
  });
  byId('assessmentClearComparison').addEventListener('click',()=>{comparisons=[];try{localStorage.removeItem(comparisonStorage);}catch{}renderComparisons();});
  byId('assessmentDownload').addEventListener('click',()=>{
    if (!report) return;
    const url=URL.createObjectURL(new Blob([JSON.stringify(report,null,2)],{type:'application/json'}));
    const a=node('a');a.href=url;a.download=`electrode_distribution_${report.job_id}.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });

  function drawCurve() {
    const box=byId('coverageCurve');box.replaceChildren();if(!report)return;
    const ns='http://www.w3.org/2000/svg', svg=document.createElementNS(ns,'svg');
    svg.setAttribute('viewBox','0 0 600 380');svg.setAttribute('aria-label','距离与面积覆盖率曲线');
    const add=(tag,attrs,text)=>{const e=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(attrs))e.setAttribute(k,String(v));if(text!=null)e.textContent=text;svg.append(e);return e;};
    const left=56,top=22,width=518,height=288,max=Math.max(2*report.target.radius_um,report.target.coverage_distance_um);
    const x=value=>left+value/max*width,y=value=>top+(1-value)*height;
    for(let i=0;i<=4;i++){const v=i/4;add('line',{x1:left,y1:y(v),x2:left+width,y2:y(v),stroke:'#e2eae4'});add('text',{x:left-8,y:y(v)+4,'text-anchor':'end'},`${Math.round(v*100)}%`);add('text',{x:x(v*max),y:top+height+22,'text-anchor':'middle'},number(v*max,0));}
    const path=(curve,key)=>curve.distance_um.map((value,i)=>`${i?'L':'M'}${x(value).toFixed(2)},${y(curve[key][i]).toFixed(2)}`).join(' ');
    const curve=report.curve;
    const lower=path(curve,'lower'),back=[...curve.distance_um].reverse().map((value,j)=>`L${x(value).toFixed(2)},${y(curve.upper[curve.upper.length-1-j]).toFixed(2)}`).join(' ');
    add('path',{d:`${lower} ${back} Z`,fill:'#61ac8d',opacity:.20});
    const others=matchingComparisons().filter(s=>s.job_id!==report.job_id).slice(-4),colors=['#8171a5','#b17b44','#5b8ca1','#9b667e'];
    others.forEach((s,i)=>add('path',{d:path(s.curve,'estimate'),fill:'none',stroke:colors[i],'stroke-width':1.5,'stroke-dasharray':'5 4'}));
    add('path',{d:path(curve,'estimate'),fill:'none',stroke:'#167866','stroke-width':2.5});
    add('line',{x1:x(report.target.coverage_distance_um),x2:x(report.target.coverage_distance_um),y1:top,y2:top+height,stroke:'#c28039','stroke-dasharray':'4 4'});
    add('text',{x:315,y:365,'text-anchor':'middle'},'到最近电极中心的距离 ℓ（μm）');
    const marker=add('circle',{r:4,fill:'#167866',visibility:'hidden'});
    svg.addEventListener('pointermove',event=>{
      const rect=svg.getBoundingClientRect(),value=Math.max(0,Math.min(max,((event.clientX-rect.left)/rect.width*600-left)/width));
      let index=0;while(index+1<curve.distance_um.length&&curve.distance_um[index+1]<=value)index++;
      marker.setAttribute('cx',x(curve.distance_um[index]));marker.setAttribute('cy',y(curve.estimate[index]));marker.setAttribute('visibility','visible');
      byId('coverageHover').textContent=`ℓ = ${number(curve.distance_um[index])} μm · 覆盖约 ${percent(curve.estimate[index])} · 区间 ${percent(curve.lower[index])}–${percent(curve.upper[index])}`;
    });svg.addEventListener('pointerleave',()=>marker.setAttribute('visibility','hidden'));
    box.append(svg);byId('coverageHover').textContent=`当前 ℓ = ${number(report.target.coverage_distance_um)} μm · 覆盖约 ${percent(report.coverage.estimate)}${others.length?'；虚线为本组其它保存布局':''}`;
    const legend=byId('coverageCurveLegend');legend.replaceChildren();
    [report,...others].forEach((item,i)=>{const label=node('span'),line=node('i');line.style.borderTop=`2px ${i?'dashed':'solid'} ${i?colors[i-1]:'#167866'}`;label.append(line,node('span',`${(item.input_name||item.job_id).split(/[\\/]/).pop()} · N=${item.electrode_count}${i?'':'（当前）'}`));legend.append(label);});
  }
  const stops=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
  function color(value){const v=Math.min(1,Math.max(0,value))*4,i=Math.min(3,Math.floor(v)),f=v-i;return stops[i].map((x,k)=>Math.round(x*(1-f)+stops[i+1][k]*f));}
  function drawHeatmap() {
    const canvas=byId('distributionHeatmap');canvas.width=canvas.height=512;const ctx=canvas.getContext('2d');
    ctx.fillStyle='#f5f8f4';ctx.fillRect(0,0,512,512);field=null;canvasBounds=null;
    const heat=report?.heatmap,t=report?.target;if(!heat){ctx.fillStyle='#6e8279';ctx.textAlign='center';ctx.font='16px sans-serif';ctx.fillText('无电极，距离未定义',256,256);byId('heatmapScale').textContent='无有限距离';return;}
    const bytes=Uint8Array.from(atob(heat.values_base64),char=>char.charCodeAt(0)),view=new DataView(bytes.buffer);
    field=new Uint16Array(heat.width*heat.height);for(let i=0;i<field.length;i++)field[i]=view.getUint16(i*2,true);
    const small=document.createElement('canvas');small.width=heat.width;small.height=heat.height;const sc=small.getContext('2d'),image=sc.createImageData(heat.width,heat.height),scale=2*t.radius_um;
    for(let row=0;row<heat.height;row++)for(let col=0;col<heat.width;col++){
      const raw=field[row*heat.width+col],index=((heat.height-1-row)*heat.width+col)*4;if(!raw)continue;
      const d=(raw-1)/65534*heat.maximum_um;image.data.set([...color(d/scale),255],index);
    }sc.putImageData(image,0,0);
    const pad=30,side=452;canvasBounds={pad,side};ctx.save();ctx.beginPath();ctx.arc(256,256,side/2,0,2*Math.PI);ctx.clip();ctx.imageSmoothingEnabled=false;ctx.drawImage(small,pad,pad,side,side);
    const xy=p=>[pad+((p[0]-t.center_um[0])/t.radius_um+1)*side/2,pad+(1-(p[1]-t.center_um[1])/t.radius_um)*side/2];
    for(const p of report.electrode_centers_um){if(Math.hypot(p[0]-t.center_um[0],p[1]-t.center_um[1])>t.radius_um)continue;const[a,b]=xy(p);ctx.beginPath();ctx.arc(a,b,3.5,0,2*Math.PI);ctx.fillStyle='#fff';ctx.fill();ctx.strokeStyle='#203c31';ctx.lineWidth=1;ctx.stroke();}
    const witness=report.maximum_uncovered.witness_um;if(witness){const[a,b]=xy(witness);ctx.strokeStyle='#ffecb8';ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(a-6,b-6);ctx.lineTo(a+6,b+6);ctx.moveTo(a-6,b+6);ctx.lineTo(a+6,b-6);ctx.stroke();}
    ctx.restore();ctx.strokeStyle='#53746a';ctx.lineWidth=1.2;ctx.beginPath();ctx.arc(256,256,side/2,0,2*Math.PI);ctx.stroke();ctx.fillStyle='#587265';ctx.font='12px sans-serif';ctx.textAlign='center';ctx.fillText(`${number((t.center_um[0]-t.radius_um)/1000)} mm`,pad,505);ctx.fillText(`${number(t.center_um[0]/1000)} mm`,256,505);ctx.fillText(`${number((t.center_um[0]+t.radius_um)/1000)} mm`,pad+side,505);
    byId('heatmapScale').textContent=`${number(scale,0)} μm（以上同色）`;byId('heatmapHover').textContent='移动鼠标查看 GDS 位置与到最近电极的距离';
  }
  byId('distributionHeatmap').addEventListener('pointermove',event=>{
    if(!field||!report||!canvasBounds)return;const rect=event.currentTarget.getBoundingClientRect(),px=(event.clientX-rect.left)/rect.width*512,py=(event.clientY-rect.top)/rect.height*512,{pad,side}=canvasBounds;
    const col=Math.floor((px-pad)/side*report.heatmap.width),row=report.heatmap.height-1-Math.floor((py-pad)/side*report.heatmap.height);
    if(col<0||row<0||col>=report.heatmap.width||row>=report.heatmap.height||Math.hypot(px-256,py-256)>side/2)return;
    const raw=field[row*report.heatmap.width+col];if(!raw)return;
    const t=report.target,x=t.center_um[0]+((px-pad)/side*2-1)*t.radius_um,y=t.center_um[1]+(1-(py-pad)/side*2)*t.radius_um,d=(raw-1)/65534*report.heatmap.maximum_um;
    byId('heatmapHover').textContent=`位置 (${number(x/1000,4)}, ${number(y/1000,4)}) mm · 距离约 ${number(d)} μm · 空间误差约 ±${number(report.numerics.distance_uncertainty_um)} μm`;
  });

  function equation(parent, tex) {
    const box=node('div',null,'formula-equation');parent.append(box);
    try { if(!window.katex)throw new Error('本地公式渲染资源尚未加载');window.katex.render(tex,box,{displayMode:true,throwOnError:true,trust:false,output:'htmlAndMathml'}); }
    catch(error){box.textContent=`公式渲染失败：${error.message}`;box.classList.add('formula-error');}
  }
  function openFormula(kind) {
    const body=byId('distributionFormulaBody');body.replaceChildren();
    const title={nearest:'最近邻间距变异系数：局部间距是否一致',maximum:'最大未覆盖距离：最远位置离电极多远',coverage:'覆盖率与曲线：指定距离能覆盖多少面积'}[kind];byId('distributionFormulaTitle').textContent=title;
    body.append(node('p','设 P 为当前已通过导出回读的全部电极中心，N 为电极数量；Ω 为评价面板独立保存的目标圆。圆外电极也计入 P，可以贡献目标圆内的覆盖；空白支撑区仍属于 Ω。'));
    equation(body,String.raw`P=\{x_1,\ldots,x_N\},\qquad\Omega=\{q\in\mathbb R^2:\|q-c\|_2\le R\}`);
    const symbols=node('table');[['xᵢ','第 i 个最终电极中心，单位 μm'],['q','目标圆内任意一个位置'],['c、R','评价目标圆的圆心与半径；独立于布线的允许放置区域'],['ℓ','评价用距离尺度，单位 μm；不是默认已验证的信号探测半径'],['Area','二维平面面积；整个目标圆面积为 πR²']].forEach(([name,description])=>{const tr=node('tr');tr.append(node('td',name),node('td',description));symbols.append(tr);});body.append(symbols);
    if(kind==='nearest'){
      body.append(node('h4','1. 每个电极到最近邻居的距离'));
      equation(body,String.raw`d_i=\min_{j\ne i}\|x_i-x_j\|_2`);
      body.append(node('h4','2. 均值、总体标准差和变异系数'));
      equation(body,String.raw`\bar d=\frac1N\sum_{i=1}^N d_i,\qquad\sigma_d=\sqrt{\frac1N\sum_{i=1}^N(d_i-\bar d)^2}`);
      equation(body,String.raw`CV_{\mathrm{NN}}=\frac{\sigma_d}{\bar d}`);
      body.append(node('p','这里用总体标准差，分母为 N；CV 没有单位。CV 越接近 0，最近邻距离越一致。但等距排成一圈也可以得到 CV=0，圆中心仍可能很空，因此必须结合 h 和覆盖率。N<2 或平均距离为 0 时，CV 未定义。'));
      body.append(node('p',report?`本布局 N=${report.electrode_count}；平均距离 ${number(report.nearest_neighbor.mean_um)} μm，标准差 ${number(report.nearest_neighbor.std_um)} μm，CV=${number(report.nearest_neighbor.cv,6)}。`:'完成评价后显示本布局代入数值。','formula-current'));
    }else{
      body.append(node('h4','1. 位置到最近电极的距离'));
      equation(body,String.raw`D(q;P)=\min_{1\le i\le N}\|q-x_i\|_2`);
      if(kind==='maximum'){
        body.append(node('h4','2. 整个目标圆内的最大值'));
        equation(body,String.raw`h(P,\Omega)=\max_{q\in\Omega}D(q;P)`);
        body.append(node('p','h 越小，目标圆里最远离电极的位置也越接近电极。它衡量最大空间空缺，以电极中心计距，不扣除电极半径。没有电极时距离视为无穷。白点是电极，热图中的 × 是当前最远采样代表点；它不一定是连续区域的精确最远点。'));
        body.append(node('h4','3. 为什么报告区间'));
        equation(body,String.raw`|D(q;P)-D(q';P)|\le\|q-q'\|_2`);
        equation(body,String.raw`\varepsilon=\frac{s}{\sqrt2},\qquad\widehat h=\max_jD(q_j;P),\qquad\widehat h\le h\le\widehat h+\varepsilon`);
        body.append(node('p','将目标圆与边长 s 的方格求交。每个格子用格心代表，圆外格心投影到圆边界。代表点与该交集内任意位置相距至多 ε；最近距离函数的变化不超过位置变化，因此给出这个空间误差区间。实现另留浮点舍入余量，显示的下端略低于采样值。提高计算精度能缩小区间。'));
        body.append(node('p',report?`当前 h ${report.maximum_uncovered.unbounded?'为无穷':`估计 ${number(report.maximum_uncovered.estimate_um)} μm，区间 ${number(report.maximum_uncovered.lower_um)}–${number(report.maximum_uncovered.upper_um)} μm`}。`:'完成评价后显示本布局代入数值。','formula-current'));
      }else{
        body.append(node('h4','2. 距离阈值内的面积比例'));
        equation(body,String.raw`C(\ell;P,\Omega)=\frac{\operatorname{Area}\{q\in\Omega:D(q;P)\le\ell\}}{\operatorname{Area}(\Omega)}`);
        equation(body,String.raw`C(\ell)=\frac1{\pi R^2}\int_\Omega\mathbf1\{D(q;P)\le\ell\}\,dq`);
        body.append(node('p','指示函数 1{条件} 在条件成立时取 1，否则取 0。C 在 0 与 1 之间；页面显示百分比。覆盖率曲线就是把 ℓ 从小到大变化：同一距离下越高，覆盖面积越大。两条曲线可能相交，不能保证某种布局在所有距离上都更好。'));
        body.append(node('h4','3. 面积加权估计与上下界'));
        equation(body,String.raw`a_j=\operatorname{Area}(\Omega\cap Q_j),\qquad\widehat C(\ell)=\frac{\sum_j a_j\mathbf1\{D(q_j;P)\le\ell\}}{\pi R^2}`);
        equation(body,String.raw`C_-(\ell)=\widehat C(\ell-\varepsilon)\le C(\ell)\le\widehat C(\ell+\varepsilon)=C_+(\ell)`);
        body.append(node('p','Qⱼ 是网格方块，aⱼ 是它与真实圆的交集面积；边界格不能和完整格算成相同面积。实现用圆积分的解析表达式计算面积，曲线阴影表示距离误差产生的上下界，并使用浮点计算。ℓ=0 时有限电极中心集的面积为 0；无电极时所有有限尺度覆盖率为 0。'));
        body.append(node('p',report?`当前 ℓ=${number(report.coverage.distance_um)} μm：覆盖率估计 ${percent(report.coverage.estimate)}，区间 ${percent(report.coverage.lower)}–${percent(report.coverage.upper)}。`:'完成评价后显示本布局代入数值。','formula-current'));
      }
    }
    body.append(node('h4','比较时保持一致'));body.append(node('p','比较项须有相同目标圆心、半径、覆盖距离及计算精度；页面会排除不匹配的保存项。同时报告总电极数和圆内电极数。此处只描述二维几何分布，不改变布局，也不证明电极数量最多或实际信号覆盖。'));
    const dialog=byId('distributionFormulaDialog');if(!dialog.open)dialog.showModal();body.scrollTop=0;
  }
  document.querySelectorAll('[data-formula]').forEach(button=>button.addEventListener('click',()=>openFormula(button.dataset.formula)));
  byId('distributionFormulaClose').addEventListener('click',()=>byId('distributionFormulaDialog').close());
  byId('distributionFormulaDialog').addEventListener('click',event=>{if(event.target===event.currentTarget){const r=event.currentTarget.getBoundingClientRect();if(event.clientX<r.left||event.clientX>r.right||event.clientY<r.top||event.clientY>r.bottom)event.currentTarget.close();}});
  window.DistributionPanel={update(nextJob,nextDetail){
    const changed=job?.id!==nextJob?.id;job=nextJob;detail=nextDetail;
    const ready=job?.status==='complete'&&!!job.result;
    byId('assessmentApply').disabled=!ready;byId('assessmentUseCenter').disabled=!Array.isArray(detail?.geometry?.electrode_region_reference_um);
    if(!ready){if(changed||report){++epoch;controller?.abort();report=null;lastKey=null;byId('distributionContent').hidden=true;}byId('distributionBadge').textContent='等待结果';status('选择一个已完成的任务，即可评价已有布局，无需重新布线。');return;}
    if(!dirty)calculate();
  }};
})();
