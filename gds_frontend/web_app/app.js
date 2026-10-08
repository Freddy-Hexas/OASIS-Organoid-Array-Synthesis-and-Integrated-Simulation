const $ = (id) => document.getElementById(id);
const state = { files: [], jobsByFile: {}, selected: null, detail: null, view: 'support', job: null, jobFullId: null, jobFullStatus: null, poll: null, pollBusy: false, submitPending: false, imageToken: 0, detailsToken: 0, stageKey: null, viewKey: null, viewer: { kind: null, svg: null, base: null, imageScale: 1, imageX: 0, imageY: 0, drag: null } };
const vectorPreview = typeof Worker !== 'undefined' && typeof OffscreenCanvas !== 'undefined' ?
  new VectorPreview($('canvasStage'), $('vectorViewport')) : null;
const VIEW_LABELS = { feasible: '紫色：原支撑内已提取的圆岛锚点域 · 矢量', graph: '自动提取的走廊、分叉与候选出口 · PNG', routing: '导出 GDS 的圆形衬底、电极与金属 · 矢量', curve: '深色：解析二次贝塞尔中心线；底图为实际导出 GDS · 矢量' };
const STAGES = [
  { view: 'support', title: '原始支撑', detail: '读取 GDS 支撑层', missing: '选择 GDS' },
  { view: 'feasible', title: '电极可放区', detail: '原结构内的锚点域', missing: '运行后生成' },
  { view: 'graph', title: '走廊与出口', detail: '自动提取的拓扑图', missing: '运行后生成' },
  { view: 'routing', title: '电极与金属', detail: '缩略图放大首条路线', missing: '运行后生成' },
  { view: 'curve', title: '解析曲线', detail: '局部放大曲线中心线', missing: '运行后生成' }
];
const supportLabel = () => `GDS ${state.detail?.selected_support_layer || '—'} 支撑层 · 原始矢量轮廓`;
const taskForFile = (file) => {
  let task = state.jobsByFile[file?.id];
  // A running server can still return its pre-update compatibility policy.
  // These explicit historical exceptions preserve results and revision;
  // changed input hashes and physically incompatible old layouts stay stale.
  if (task?.stale_reason === 'pad_layout_algorithm_changed' &&
      ['center_preferred_pad_routing_v12_electrode_region',
       'center_preferred_pad_routing_v13_attachment_cells'].includes(task.solver_revision)) {
    task = { ...task, stale_reason: null,
      parameter_note: task.solver_revision === 'center_preferred_pad_routing_v13_attachment_cells' ?
        '历史布局采用限时搜索；新任务取消默认求解时间限制。原结果保留供查看，尚未按不限时策略重跑。' :
        '历史布局采用旧位置搜索；新任务使用圆岛附着自适应搜索。原结果保留供查看，尚未按新策略重跑。' };
  }
  return task?.stale_reason || (task?.input_sha256 && file?.sha256 && task.input_sha256 !== file.sha256) ? null : task;
};
const jobMatches = () => state.job?.file_id === state.detail?.id &&
  (state.job?.support_layer || '10/0') === state.detail?.selected_support_layer &&
  (!state.job?.result?.geometry?.sha256 || state.job.result.geometry.sha256 === state.detail?.sha256);
const historicalResultNote = () => {
  const task = taskForFile(state.selected);
  return task && state.job && task.id === state.job.id ? task.parameter_note || '' : '';
};
const currentSummary = () => jobMatches() && state.job?.result ? state.job.result : state.detail?.reference?.summary;
const padRun = (result) => result?.method === 'four_side_pads';

function toast(message, error = false) {
  const node = $('toast');
  node.textContent = message;
  node.classList.toggle('error', error);
  node.classList.add('show');
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => node.classList.remove('show'), 3600);
}

async function getJSON(url, options) {
  const response = await fetch(url, { cache: 'no-store', signal: AbortSignal.timeout(15000), ...options });
  let data;
  try { data = await response.json(); } catch { data = {}; }
  if (!response.ok) throw new Error(data.error || `HTTP ${response.status}`);
  return data;
}

const fmt = (value) => value == null ? '—' : Number(value).toLocaleString('zh-CN');
const shortBytes = (value) => value < 1e6 ? `${Math.round(value / 1024)} KB` : `${(value / 1e6).toFixed(2)} MB`;
function terminalBreakdown(result) {
  if (result?.method !== 'attached_islands' || !Array.isArray(result.routing?.routes)) return null;
  const counts = { collector_interface: 0, open_tip: 0, unknown: 0 };
  for (const route of result.routing.routes) {
    const kind = route.terminal_kind || (Number.isFinite(result.geometry?.triangle_count) &&
      Number.isFinite(route.outlet_node) ?
      (route.outlet_node >= result.geometry.triangle_count ? 'collector_interface' : 'open_tip') : 'unknown');
    counts[kind in counts ? kind : 'unknown'] += 1;
  }
  return counts;
}
function terminalSummary(result) {
  const counts = terminalBreakdown(result);
  if (!counts) return '';
  const parts = [`${fmt(counts.collector_interface)} 条止于宽区接口`, `${fmt(counts.open_tip)} 条止于开放末端`];
  if (counts.unknown) parts.push(`${fmt(counts.unknown)} 条终点未分类`);
  return parts.join('、');
}

function setStateBadge(text, style = '') {
  const node = $('fileState');
  node.textContent = text;
  node.className = `state-pill ${style}`;
}

function create(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text != null) node.textContent = text;
  return node;
}

function renderFiles() {
  const box = $('fileList'); box.replaceChildren();
  const query = $('searchInput').value.trim().toLowerCase();
  const visible = state.files.filter(f => f.routable_input &&
    `${f.name} ${f.relative_path}`.toLowerCase().includes(query));
  if (!visible.length) { box.append(create('p', 'muted-copy', '没有匹配的 GDS 文件。')); return; }
  visible.forEach(file => {
    const button = create('button', 'file-row' + (state.selected?.id === file.id ? ' active' : ''));
    button.type = 'button'; button.setAttribute('role', 'option');
    button.setAttribute('aria-selected', state.selected?.id === file.id ? 'true' : 'false');
    button.append(create('span', 'file-icon', 'GDS'));
    const copy = create('span', 'file-copy'); copy.append(create('span', 'file-name', file.name));
    copy.append(create('span', 'file-sub', file.relative_path === file.name ? `${shortBytes(file.size_bytes)} · 根目录` : file.relative_path));
    const task = taskForFile(file);
    const stale = !!state.jobsByFile[file.id] && !task;
    const completed = task?.method === 'four_side_pads' ?
      `${fmt(task.lower_bound ?? 0)} 个电极接通 Pad` :
      task?.lower_bound == null ? '已分析 · 无候选终端路线' : `试布线 ${task.lower_bound} 个 · 未接 Pad`;
    const status = stale ?
      (state.jobsByFile[file.id].stale_reason === 'pad_layout_algorithm_changed' ? 'Pad 算法已更新 · 请重跑' : 'GDS 已变化 · 请重跑') :
      task ? ({ complete: completed, error: '运行失败', running: '运行中', queued: '排队中', interrupted: '已中断' })[task.status] || task.status : '待运行';
    copy.append(create('span', `task-status ${task?.status || ''}`, status));
    if (task?.status === 'complete' && task.rules?.minimum_center_spacing_um != null)
      copy.append(create('span', 'task-rule', `中心距 ${task.rules.minimum_center_spacing_um === 0 ? '未设置' : `${fmt(task.rules.minimum_center_spacing_um / 1000)} mm`}`));
    if (task?.parameter_note) copy.append(create('span', 'task-rule', '历史布局 · 新策略尚未重跑'));
    button.append(copy);
    if (task?.status === 'complete' && task.lower_bound != null) {
      const count = create('span', 'file-electrode-count');
      count.append(create('strong', '', String(task.lower_bound)), create('small', '', '电极'));
      button.append(count);
    }
    button.addEventListener('click', () => selectFile(file.id));
    box.append(button);
  });
}

function fileStatusSignature() {
  return state.files.map(file => {
    const task = taskForFile(file);
    return `${file.id}:${task?.id || ''}:${task?.status || ''}:${task?.lower_bound ?? ''}:${task?.stale_reason || ''}:${task?.parameter_note || ''}`;
  }).join('|');
}

async function refreshFiles() {
  $('refreshButton').disabled = true;
  try {
    const data = await getJSON('/api/files');
    state.files = data.files.filter(f => f.routable_input).sort((a, b) =>
      a.name.localeCompare(b.name, 'zh-CN'));
    state.jobsByFile = data.tasks || {};
    $('rootPath').textContent = data.root;
    $('fileCount').textContent = String(state.files.length).padStart(2, '0');
    $('lastScan').textContent = new Date().toLocaleTimeString('zh-CN', { hour: '2-digit', minute: '2-digit' });
    renderFiles();
    const requestedId = new URLSearchParams(location.search).get('file');
    const requested = state.files.find(f => f.id === requestedId);
    const preferred = state.selected && state.files.find(f => f.id === state.selected.id);
    const savedFile = state.files.find(f => f.id === localStorage.getItem('gds-workbench-file'));
    const example = state.files.find(f => f.name === 'C_open_petal_mesh.gds');
    const first = requested || preferred || savedFile || example || state.files[0];
    if (first) await selectFile(first.id, first.id === preferred?.id ? state.detail?.selected_support_layer : null);
    if (!state.poll) state.poll = setInterval(pollQueue, 1500);
    pollQueue();
    if (!data.files.length) toast('固定目录下没有 GDS 文件。', true);
  } catch (error) { toast(`文件扫描失败：${error.message}`, true); }
  finally { $('refreshButton').disabled = false; }
}

function resetZoom() {
  const v = state.viewer;
  if (v.kind === 'vector') vectorPreview.reset();
  if (v.kind === 'svg' && v.svg && v.base) v.svg.setAttribute('viewBox', v.base.join(' '));
  if (v.kind === 'image') {
    v.imageScale = 1; v.imageX = 0; v.imageY = 0;
    $('mainImage').style.transform = 'translate(0px, 0px) scale(1)';
  }
}

function cancelDrag(pointerId = null) {
  const drag = state.viewer.drag;
  if (pointerId != null && drag?.pointerId !== pointerId) return;
  state.viewer.drag = null;
  const stage = $('canvasStage');
  stage.classList.remove('dragging');
  if (drag && stage.hasPointerCapture(drag.pointerId)) stage.releasePointerCapture(drag.pointerId);
}

function svgPoint(clientX, clientY) {
  const svg = state.viewer.svg;
  const point = svg.createSVGPoint(); point.x = clientX; point.y = clientY;
  return point.matrixTransform(svg.getScreenCTM().inverse());
}

function zoom(factor, clientX, clientY) {
  const v = state.viewer;
  if (!v.kind) return;
  if (v.kind === 'vector') { vectorPreview.zoom(factor, clientX, clientY); return; }
  const rect = $('canvasStage').getBoundingClientRect();
  clientX ??= rect.left + rect.width / 2; clientY ??= rect.top + rect.height / 2;
  if (v.kind === 'svg') {
    const svg = v.svg, old = svg.viewBox.baseVal;
    const base = v.base;
    const newWidth = Math.max(base[2] / 500, Math.min(base[2] * 2, old.width / factor));
    const ratio = newWidth / old.width;
    const point = svgPoint(clientX, clientY);
    svg.setAttribute('viewBox', `${point.x - (point.x - old.x) * ratio} ${point.y - (point.y - old.y) * ratio} ${newWidth} ${old.height * ratio}`);
  } else {
    const next = Math.max(1, Math.min(40, v.imageScale * factor));
    const actual = next / v.imageScale;
    const dx = clientX - (rect.left + rect.width / 2), dy = clientY - (rect.top + rect.height / 2);
    v.imageX = dx - (dx - v.imageX) * actual;
    v.imageY = dy - (dy - v.imageY) * actual;
    v.imageScale = next;
    $('mainImage').style.transform = `translate(${v.imageX}px, ${v.imageY}px) scale(${v.imageScale})`;
  }
}

async function setView(url, alt, vector) {
  const token = ++state.imageToken;
  const img = $('mainImage'), viewport = $('vectorViewport');
  cancelDrag();
  vectorPreview?.cancel();
  state.viewAbort?.abort(); state.viewAbort = new AbortController();
  state.viewer = { kind: null, svg: null, base: null, imageScale: 1, imageX: 0, imageY: 0, drag: null };
  img.hidden = true; viewport.hidden = true; viewport.replaceChildren();
  $('imageLoading').hidden = false; $('canvasEmpty').hidden = true;
  try {
    if (vector) {
      if (vectorPreview) {
        viewport.hidden = false;
        await vectorPreview.load(url.replace('/api/vector/', '/api/scene/'), alt);
        if (token !== state.imageToken) return;
        state.viewer = { kind: 'vector', drag: null };
        updateElectrodeRegion();
        $('imageLoading').hidden = true;
        renderThumbnails();
        return;
      }
      const response = await fetch(url, { signal: state.viewAbort.signal });
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const documentSVG = new DOMParser().parseFromString(await response.text(), 'image/svg+xml');
      if (documentSVG.querySelector('parsererror')) throw new Error('SVG 解析失败');
      const svg = documentSVG.documentElement;
      if (svg.localName !== 'svg') throw new Error('不是 SVG');
      if (token !== state.imageToken) return;
      viewport.append(document.importNode(svg, true));
      const shown = viewport.querySelector('svg');
      state.viewer = { kind: 'svg', svg: shown, base: shown.getAttribute('viewBox').split(/\s+/).map(Number), drag: null };
      viewport.hidden = false;
    } else {
      const loaded = new Image();
      await new Promise((resolve, reject) => {
        loaded.onload = resolve;
        loaded.onerror = () => reject(new Error('PNG 读取失败'));
        loaded.src = url;
      });
      if (token !== state.imageToken) return;
      img.src = loaded.src;
      img.alt = alt; img.hidden = false;
      state.viewer.kind = 'image'; resetZoom();
    }
    $('imageLoading').hidden = true;
    updateElectrodeRegion();
    renderThumbnails();
  } catch (error) {
    if (token !== state.imageToken) return;
    if (error.name === 'AbortError') return;
    state.viewKey = null;
    $('imageLoading').hidden = true; $('canvasEmpty').hidden = false;
    toast(`视图读取失败：${error.message}`, true);
  }
}

function viewURL(view = state.view) {
  const detail = state.detail;
  if (!detail) return null;
  if (view === 'support') return { url: `/api/vector/input/${detail.id}?spec=${encodeURIComponent(detail.selected_support_layer)}&v=${detail.sha256}`, vector: true };
  if (view === 'feasible') {
    if (jobMatches() && state.job?.artifacts?.includes('regions.json.gz')) return { url: `/api/vector/feasible/${state.job.id}`, vector: true };
    return null;
  }
  if (view === 'graph') {
    if (jobMatches() && state.job?.artifacts?.includes('graph.png')) return { url: `/api/artifact/${state.job.id}/graph.png`, vector: false };
    return detail.reference?.graph_image ? { url: detail.reference.graph_image, vector: false } : null;
  }
  if (view === 'routing') {
    if (jobMatches() && state.job?.artifacts?.includes('routing.gds')) return { url: `/api/vector/job/${state.job.id}`, vector: true };
    return detail.reference?.routing_image ? { url: detail.reference.routing_image, vector: false } : null;
  }
  if (view === 'curve' && jobMatches() && (state.job?.result?.routing?.preview_metrics?.has_curves || state.job?.result?.routing?.routes?.some(r => r.curve?.curve_segments?.length)))
    return { url: `/api/vector/curve/${state.job.id}`, vector: true };
  return null;
}

function renderStages() {
  const job = jobMatches() ? state.job : null;
  const cornerCount = job?.result?.routing?.preview_metrics?.curved_corner_count ??
    job?.result?.routing?.routes?.reduce((sum, route) => sum + (route.curve?.curved_corner_count || 0), 0) ?? 0;
  const curved = cornerCount > 0;
  const stats = job?.result || state.detail?.reference?.summary;
  const geometry = stats?.geometry || state.detail?.geometry;
  const key = [state.detail?.id, state.detail?.selected_support_layer, job?.id,
    (job?.artifacts || []).join(','), curved, state.detail?.reference?.input_name].join('|');
  const box = $('stageCards');
  if (state.stageKey !== key) {
    state.stageKey = key;
    box.replaceChildren();
    STAGES.forEach((stage, index) => {
      const source = viewURL(stage.view);
      const card = create('button', 'stage-card');
      card.type = 'button'; card.dataset.view = stage.view; card.disabled = !source;
      const top = create('span', 'stage-card-top');
      const readyLabel = stage.view === 'support' ?
        (geometry?.support_holes == null ? '原始 GDS' : `${fmt(geometry.support_holes)} 个孔洞`) :
        stage.view === 'feasible' && stats?.process_geometry?.attachment_anchor_region_area_um2 != null ?
          `${fmt(Math.round(stats.process_geometry.attachment_anchor_region_area_um2))} μm²` :
        stage.view === 'graph' && (geometry?.geometry_candidate_terminals ?? geometry?.collector_interface_windows) != null ?
          `${fmt(geometry.geometry_candidate_terminals ?? geometry.collector_interface_windows)} 个出口` :
        stage.view === 'routing' && job?.result?.routing?.gds_roundtrip_audit?.passed ? `${job.result.routing.retained_routes} 个电极` :
        stage.view === 'curve' && curved ? `${cornerCount} 个圆角` :
        source?.url.includes('/api/reference/') ? '历史参考' : '已生成';
      top.append(create('span', 'stage-number', String(index + 1).padStart(2, '0')),
        create('span', `stage-status ${source ? 'ready' : ''}`,
          source ? readyLabel :
            job?.status === 'complete' && ['routing', 'curve'].includes(stage.view) ? '无完整路线' : stage.missing));
      const thumbnail = create('span', `stage-thumbnail ${source ? '' : 'unavailable'}`);
      if (source) {
        const image = create(source.vector && vectorPreview ? 'canvas' : 'img');
        const preview = stage.view === 'routing' && job?.artifacts?.includes('routing.gds') ?
          `/api/vector/job/${job.id}?highlight=metal&focus=first_route` :
          stage.view === 'curve' && curved ? `/api/vector/curve/${job.id}?focus=first_route` : source.url;
        image.alt = ''; image.loading = 'lazy'; image.decoding = 'async';
        const loading = create('span', 'stage-thumb-loading', '载入预览…');
        image.addEventListener('load', () => thumbnail.classList.add('loaded'));
        image.addEventListener('error', () => { loading.textContent = '预览暂不可用'; thumbnail.classList.add('failed'); });
        thumbnail.append(loading, image);
        if (source.vector && vectorPreview) {
          image.dataset.scene = source.url.replace('/api/vector/', '/api/scene/');
          image.dataset.focus = ['routing', 'curve'].includes(stage.view) ? 'true' : '';
          image.dataset.highlight = stage.view === 'routing' ? 'true' : '';
        } else image.src = preview;
      } else {
        thumbnail.append(create('span', 'stage-placeholder', index === 0 ? 'GDS' : '—'));
      }
      const copy = create('span', 'stage-copy');
      copy.append(create('strong', '', stage.title), create('small', '', stage.detail));
      card.append(top, thumbnail, copy);
      if (source) card.addEventListener('click', () => { state.view = stage.view; renderView(); });
      box.append(card);
    });
  }
  box.querySelectorAll('.stage-card').forEach(card => {
    const active = card.dataset.view === state.view;
    card.classList.toggle('active', active);
    card.setAttribute('aria-pressed', String(active));
  });
}

function renderThumbnails() {
  if (!vectorPreview) return;
  $('stageCards').querySelectorAll('canvas[data-scene]').forEach(canvas => {
    if (canvas.closest('.stage-thumbnail').classList.contains('loaded') ||
        canvas.dataset.queuedToken === String(vectorPreview.token)) return;
    canvas.dataset.queuedToken = String(vectorPreview.token);
    vectorPreview.thumbnail(canvas.dataset.scene, canvas,
      { focus: !!canvas.dataset.focus, highlight: !!canvas.dataset.highlight });
  });
}

function renderView() {
  renderStages();
  $('canvasTitle').textContent = STAGES.find(stage => stage.view === state.view)?.title || '几何预览';
  const label = state.view === 'support' ? supportLabel() :
    state.view === 'routing' && padRun(currentSummary()) ? '真实 GDS：电极、承载桥、四边 Pad 和完整金属 · 矢量' : VIEW_LABELS[state.view];
  $('viewDescription').textContent = label;
  $('canvasSource').textContent = state.view === 'support' ? `GDS VECTOR · ${state.detail?.selected_support_layer || '—'}` : state.view === 'feasible' ? 'VECTOR REGION · CURRENT RUN' : state.view === 'curve' ? 'GDS + ANALYTIC CURVE · CURRENT RUN' : jobMatches() && state.job?.artifacts?.includes(state.view === 'graph' ? 'graph.png' : 'routing.gds') ? (state.view === 'routing' ? 'GDS VECTOR · CURRENT RUN' : 'CURRENT RUN · PNG') : state.detail?.reference ? 'REFERENCE · PNG' : '—';
  const url = viewURL();
  if (url) {
    const key = JSON.stringify([url.url, url.vector, state.detail.sha256]);
    if (key !== state.viewKey) {
      state.viewKey = key;
      setView(url.url, `${state.detail.name} · ${label}`, url.vector);
    }
  }
  else {
    state.viewKey = null;
    cancelDrag();
    vectorPreview?.cancel(); state.viewAbort?.abort();
    ++state.imageToken;
    state.viewer.kind = null; $('mainImage').hidden = true; $('vectorViewport').hidden = true; $('imageLoading').hidden = true; $('canvasEmpty').hidden = false;
    $('canvasEmpty').querySelector('strong').textContent = state.view === 'support' ? '等待结构文件' : '此视图尚未生成';
    $('canvasEmpty').querySelector('span:last-child').textContent = state.view === 'support' ? '选择 GDS 后展示真实支撑层' :
      state.view === 'feasible' ? '运行几何分析后显示电极可放区' : state.view === 'graph' ? '运行几何分析后显示自动提取图' : '运行试布线后显示金属路径';
  }
}

function renderDetails() {
  const d = state.detail;
  if (!d) return;
  $('selectedTitle').textContent = d.name;
  $('selectedSubtitle').textContent = d.relative_path;
  $('pathLabel').textContent = d.absolute_path;
  $('sizeLabel').textContent = shortBytes(d.size_bytes);
  $('hashLabel').textContent = d.sha256;
  $('copyHash').disabled = false;
  const summary = jobMatches() && state.job && ['queued', 'running'].includes(state.job.status) ? null : currentSummary();
  const geom = summary?.geometry || d.geometry;
  $('holeMetric').textContent = fmt(geom?.support_holes);
  $('componentMetric').textContent = fmt(geom?.support_components);
  $('edgeMetric').textContent = fmt(geom?.graph_edges);
  $('outletMetric').textContent = fmt(geom?.geometry_candidate_terminals ?? geom?.collector_interface_windows);
  $('detailStatus').textContent = summary?.geometry ?
    (summary.method === 'attached_islands' || padRun(summary) ? '矢量拓扑提图' : `提图 ${geom.selected_pitch_um} μm`) : '读取原始几何';
  const layers = Object.entries(d.layers || {});
  $('layerCount').textContent = `${layers.length} 类`;
  $('layerList').replaceChildren();
  layers.forEach(([layer, count]) => $('layerList').append(create('span', 'layer-chip' + (layer === d.selected_support_layer ? ' selected' : ''), `${layer} · ${fmt(count)}`)));
  const selector = $('supportLayer'); selector.replaceChildren();
  layers.forEach(([layer, count]) => {
    const option = create('option', '', `${layer} · ${fmt(count)} 多边形`);
    option.value = layer; selector.append(option);
  });
  selector.value = d.selected_support_layer; selector.disabled = false;
  $('layerReason').textContent = d.selected_support_layer === d.suggested_support_layer ?
    `自动建议 ${d.suggested_support_layer}：${d.layer_suggestion_reason}。` :
    `手动选择 ${d.selected_support_layer}；自动建议为 ${d.suggested_support_layer}。请确认该层确实是支撑材料。`;
  const ownTask = taskForFile(d);
  if (!d.routable_input) setStateBadge('历史演示 · 只读', 'warn');
  else if (ownTask?.status === 'queued') setStateBadge('已入队 · 等待运行');
  else if (ownTask?.status === 'running') setStateBadge('正在运行 · ' + (ownTask.stage || '计算中'));
  else if (!summary) setStateBadge('待自动提图');
  else if (summary.routing?.status === 'no_sampled_legal_electrode_candidate') setStateBadge('无可用采样电极', 'warn');
  else if (summary.routing?.status === 'no_sampled_route_to_selected_outlet') setStateBadge('未找到合法路径', 'warn');
  else if (padRun(summary) && summary.routing?.pad_connection_verified) setStateBadge(`${fmt(summary.routing.retained_routes)} 个 Pad 已接通`, 'good');
  else if ((geom?.geometry_candidate_terminals ?? geom?.collector_interface_windows) > 0) setStateBadge('检测到候选出口', 'good');
  else setStateBadge('尚无可确认出口', 'warn');
  const ownActive = ownTask && ['queued', 'running'].includes(ownTask.status);
  $('runButton').disabled = !d.routable_input || !!ownActive || state.submitPending;
  $('runButton').querySelector('span').textContent = ownActive ?
    (ownTask.status === 'queued' ? '此 case 已排队' : '此 case 正在运行') : '开始分析与试布线';
  const fullPad = $('routingMode').value === 'four_side_pads';
  $('runBadge').textContent = fullPad ? '四边 Pad 联合布线' : '几何终端试布线';
  $('padSettingsPanel').hidden = !fullPad;
  $('collectorClearanceField').hidden = fullPad;
  if (fullPad) {
    if ($('outletMode').value !== 'external_boundary') state.lastGeometryOutlet = $('outletMode').value;
    $('outletMode').value = 'external_boundary';
    $('outletMode').disabled = true;
    $('outletHelp').textContent = '从实际 GDS 的方向外包络提取出口：方形的四条边、圆形的外周均可形成窗口。整条承载桥必须避开内缩包络的内部区域，只允许一次接出，每条网络连接独立 Pad。';
  } else {
    $('outletMode').disabled = false;
    if ($('outletMode').value === 'external_boundary') $('outletMode').value = state.lastGeometryOutlet || 'auto_geometry';
  }
  if (!fullPad) {
    if (summary && geom?.collector_interface_windows === 0 && $('outletMode').value === 'collector') $('outletHelp').textContent = '当前 GDS 无宽汇集区；可选择实验性开放端点模式。';
    else if ($('outletMode').value === 'open_tips') $('outletHelp').textContent = '把结构开放末端作候选出口；不等于已有焊盘或扇出。';
    else if ($('outletMode').value === 'auto_geometry') $('outletHelp').textContent = '统一提取宽汇集区接口与开放末端；它们都是几何候选出口，不等于焊盘。';
    else $('outletHelp').textContent = '自动宽汇集区接口为候选出口；不等于焊盘。';
  }
  renderAnalysis();
}

function renderAnalysis() {
  const d = state.detail;
  $('analysisEmpty').hidden = !!d;
  $('analysisContent').hidden = !d;
  if (!d) { $('analysisBadge').textContent = '等待文件'; return; }
  const activeJob = jobMatches() && state.job && ['queued', 'running'].includes(state.job.status);
  const summary = activeJob ? null : currentSummary();
  const g = summary?.geometry || d.geometry;
  const p = summary?.process_geometry;
  const r = summary?.routing;
  const cap = summary?.capacity_interval;
  const isPad = padRun(summary) || (activeJob && state.job.method === 'four_side_pads');
  const isVector = summary?.method === 'attached_islands' || isPad;
  const layerKeys = Object.keys(d.layers || {});
  $('analysisBadge').textContent = activeJob ? '当前运行中' : historicalResultNote() ? '历史布局 · 旧位置策略' : jobMatches() && state.job?.result ? '当前运行结果' : summary ? '历史基线结果' : '原始 GDS 阅读';
  $('supportExplanation').textContent = `读取实际 GDS 的 ${layerKeys.length} 类图层；当前按 ${d.selected_support_layer} 解释为支撑材料。支撑面积 ${fmt(Math.round(g.support_area_um2 || 0))} μm²，包含 ${fmt(g.support_components)} 个连通块和 ${fmt(g.support_holes)} 个孔洞。${d.selected_support_layer === d.suggested_support_layer ? '这是自动建议层，未知来源仍需确认材料语义。' : '这是手动选层。'}`;
  $('topologyExplanation').textContent = activeJob ? '正在按当前工艺参数重新提取拓扑与走廊。' : g.full_graph_matches_vector_topology ?
    (isVector ?
      `原支撑矢量内缩后做约束三角剖分；三角形对偶图的连通数和环路秩与矢量走线域一致。得到 ${fmt(g.graph_nodes)} 个节点、${fmt(g.graph_edges)} 条连接；不依赖结构名称或生成器中心线。` :
      `自动骨架图使用 ${g.selected_pitch_um} μm 栅格，环路秩 ${g.full_graph_cycle_rank} 与矢量孔洞数一致；得到 ${fmt(g.graph_nodes)} 个节点、${fmt(g.graph_edges)} 条走廊。${g.raster_repair ? `栅格化产生假拓扑，已用 ${g.raster_repair.radius_um} μm 单像素内缩修复，并确认骨架点位于原始矢量支撑内。` : ''}拓扑一致仍需逐走廊矢量检查。`) :
    '尚未对本文件运行自动提图；选择工艺参数并启动分析后生成。';
  $('placementExplanation').textContent = activeJob && isPad ?
    '正在从原始 GDS 的方向外包络提取真实外端，联合分配圆岛电极、内部路径和出口。电极中心留在原支撑，圆形承载岛可局部扩充；内部导线必须留在原结构内。' : isVector ?
    (isPad ? `在原支撑的完整走线域中选圆岛锚点。采样 ${fmt(g.external_boundary_terminals)} 个外端候选门户，搜索前核查原结构接点与整条桥；${r?.outer_exit_policy?.maximum_launch_depth_from_envelope_um != null ? `从各个方向的原结构外包络内 ${fmt(r.outer_exit_policy.maximum_launch_depth_from_envelope_um)} μm 的工艺带接出，桥不得进入内缩包络核心。` : r?.outer_exit_policy ? `本历史结果采用允许接出半径至少 ${fmt(r.outer_exit_policy.minimum_launch_radius_um / 1000)} mm 的旧径向规则。` : ''}联合分配内部路径与出口，再安排连续 Pad；内部路线沿原结构走到外端。` :
    `沿原支撑内、扣除宽汇集区的走线域选圆岛锚点；每个 30 μm 电极外加承载圆岛并核查单片连接、孔洞和相邻结构距离。识别 ${fmt(g.geometry_candidate_terminals ?? g.collector_interface_windows)} 个几何终端；它们不是焊盘。`) : p ?
    `在完整支撑上按 ${p.rules.electrode_diameter_um} μm 电极直径和 ${p.rules.margin_um} μm 余量求中心可行区。识别 ${fmt(g.collector_interface_windows)} 个候选宽汇集区接口；这些接口不是焊盘。` :
    '候选电极必须同时满足实体圆盘包含、支撑边缘余量与引线连通。宽汇集区接口仅作候选出口，不能视为已定义焊盘。';
  const centerStats = r?.center_preference?.selected;
  const allowedRadius = (activeJob ? state.job?.rules : r?.rules)?.electrode_region_radius_um;
  if (allowedRadius != null) $('placementExplanation').textContent +=
    ` 电极中心另须在结构几何圆心的 ${fmt(allowedRadius / 1000)} mm 半径圆内；走线和 Pad 可在该圆外。`;
  const compaction = r?.center_compaction;
  const attachmentSearchText = compaction?.placement_search_revision ?
    '向内移动同时检查圆岛附着，并自适应细分、覆盖不同位置区域；单点失败不会排除整段结构。仍未完成细分的区域保留在报告中，位置最优性尚未证明。' : '';
  const timePolicyText = compaction?.wall_clock_limit_enabled === false ?
    '本次位置优化不设时间上限，按候选、细分精度和轮次完成；已用时间只作统计。' : '';
  const augmentation = r?.port_augmentation;
  const joint = r?.joint_path_flow;
  const jointText = joint?.individually_certified_path_proposals != null ? `联合骨干路径 ${fmt(joint.guide_flow_count)} 条，路径与电极接入位置的几何合格候选 ${fmt(joint.individually_certified_path_proposals)} 个${joint.suffix_contact_proposals != null ? `，包含 ${fmt(joint.suffix_contact_proposals)} 个沿同一路径重新选择的接入位置` : ''}；真实金属冲突联合选择保留 ${fmt(joint.certified_internal_count)} 条内部路径，初次接通 Pad ${fmt(joint.certified_complete_count)} 条。以上是提案阶段数量，最终电极数以导出审计的完整网络为准。` : '';
  const sizing = r?.pad_sizing;
  const sizingText = r?.pad_settings ? `实际 Pad 尺寸 ${fmt(r.pad_settings.pad_width_um / 1000)} × ${fmt(r.pad_settings.pad_length_um / 1000)} mm，中心节距 ${fmt(r.pad_settings.pad_pitch_um / 1000)} mm。${sizing?.pad_dimensions_reduced ? '已根据数量自动缩小，未低于设定下限。' : ''}${sizing?.frame_expanded ? '外框已按本次数量自动扩大。' : ''}` : '';
  const ports = r?.outer_port_coverage;
  const portText = ports ? `识别 ${fmt(ports.geometric_window_count)} 个几何外端窗口，已用 ${fmt(ports.used_window_count)} 个；一个窗口可容纳多轨。${augmentation ? `动态补充内部路线并共同重排外部 Pad，将完整网络从 ${fmt(augmentation.initial_connected_count)} 条增至 ${fmt(augmentation.final_connected_count)} 条。` : ''}${ports.unused_window_count ? `剩余 ${fmt(ports.unused_window_count)} 个窗口未被本次搜索使用，这不是不可布线证明。` : ''}` : '';
  $('overviewPortCoverage').textContent = ports ? `${fmt(ports.used_window_count)} / ${fmt(ports.geometric_window_count)}` : '待计算';
  const centerText = centerStats?.mean_radius_um != null ? `先在候选库中增加接通数，再向原支撑最小外接圆的中心优化；本次电极平均距中心 ${fmt(Math.round(centerStats.mean_radius_um))} μm，最远 ${fmt(Math.round(centerStats.maximum_radius_um))} μm，位于支撑最大半径一半以内 ${fmt(centerStats.within_half_source_radius)} 个。${compaction ? `残余空间重布线保留 ${fmt(compaction.before.count)} 条网络，接受 ${fmt(compaction.accepted_update_count ?? compaction.accepted_updates?.length ?? 0)} 次向内移动，平均半径从 ${fmt(Math.round(compaction.before.mean_radius_um))} 降至 ${fmt(Math.round(compaction.after.mean_radius_um))} μm；Pad、桥和外端出口固定。` : ''}` : '';
  if (activeJob) $('routingExplanation').textContent = '正在试布线；完成导出回读后再给出本次下界。';
  else if (!r) $('routingExplanation').textContent = '尚未试布线。运行后以实际金属多边形和回读 GDS 给出构造性下界。';
  else if (isPad && r.pad_connection_verified) $('routingExplanation').textContent =
    `已将 ${fmt(r.retained_routes)} 个电极分别接到四边实体 Pad，并在导出 GDS 回读后核查每条网络的连续性、承载桥、接触、边缘余量和网络间距。${centerText}${r.pad_bank_occupancy ? `四边已接通 Pad：上 ${fmt(r.pad_bank_occupancy.top.connected_pad_count)}、右 ${fmt(r.pad_bank_occupancy.right.connected_pad_count)}、下 ${fmt(r.pad_bank_occupancy.bottom.connected_pad_count)}、左 ${fmt(r.pad_bank_occupancy.left.connected_pad_count)}；每边连续无空槽，相邻 Pad 金属边缘相距 ${fmt(r.pad_settings.pad_pitch_um - r.pad_settings.pad_width_um)} μm。` : ''}本次外框边长 ${fmt(r.frame?.square_side_um / 1000)} mm、每边可用槽位 ${fmt(r.frame?.pad_slots_per_side)} 个；槽位数会随外框增大，只导出已接通的 Pad。严格整数审计确认的下界 L 为 ${fmt(cap?.lower_bound)}，连续容量上界 U 为 ${fmt(cap?.upper_bound)}。${cap?.declared_geometric_model_optimality_proven ? '两界相等，已证明当前几何模型内最多；制造 DRC 仍需单独认证。' : '两界尚未闭合，不能声称最多。'}`;
  else if (isPad && r.status === 'no_finite_complete_pad_route') $('routingExplanation').textContent =
    '本次有限候选没有找到电极到实体 Pad 的完整网络。可查看报告中的门户逃逸、桥和 Pad 接入排除原因；0 不是连续几何不可行证明。';
  else if (r.status === 'checked_geometry_lower_bound') $('routingExplanation').textContent =
    `找到 ${fmt(r.retained_routes)} 条独立金属路径，回读检查${r.gds_roundtrip_audit?.passed ? '通过' : '未完成'}，检查范围是电极到几何候选终端。${terminalSummary(summary) ? `${terminalSummary(summary)}。` : ''}最小支撑边缘净距 ${r.gds_roundtrip_audit?.minimum_metal_to_support_boundary_um?.toFixed(3) ?? '—'} μm。${r.candidate_library?.centerline_model ? '内部转角已改为与相邻直段相切的曲线；解析曲线可在专用视图查看。' : ''}${summary?.method === 'attached_islands' ? `当前锚点域的连续几何上界 ${fmt(cap?.upper_bound)}${r.center_distance_upper_bound?.angular_upper_bound != null ? `，其中角度界为 ${fmt(r.center_distance_upper_bound.angular_upper_bound)}` : ''}；衬底圆岛已并入输出支撑层，导线仍限定在原支撑。尚未认证真实焊盘。` : '连续容量上界仍未证明。'}${r.outlet_mode === 'open_tips' ? '开放端点为实验性出口。' : ''}`;
  else if (r.status === 'no_sampled_legal_electrode_candidate') $('routingExplanation').textContent = `当前有限采样位置中，没有一个能完整放下 ${p?.rules?.electrode_diameter_um ?? '—'} μm 电极并保留 ${p?.rules?.margin_um ?? '—'} μm 边缘余量。紫色“电极可放区”仍可能存在；这不是连续几何不可行证明。`;
  else if (r.status === 'no_sampled_route_to_selected_outlet') $('routingExplanation').textContent = '采样单轨图中未找到可达路线；这不证明连续几何不可布线。';
  else $('routingExplanation').textContent = '本文件没有所选模式的几何候选出口，因此只输出结构分析，不给布线数量。';
  $('electrodeAreaLabel').textContent = isVector ? '圆岛锚点域' : '电极中心可行区';
  if (portText || jointText || sizingText || attachmentSearchText || timePolicyText) $('routingExplanation').textContent = portText + jointText + sizingText + attachmentSearchText + timePolicyText + $('routingExplanation').textContent;
  $('legalCorridorsLabel').textContent = isPad ? '完整 Pad 候选' : isVector ? '已认证路径列' : '合法走廊';
  $('electrodeArea').textContent = isVector && p ? `${fmt(Math.round(p.attachment_anchor_region_area_um2))} μm²（锚点）` : p ? `${fmt(Math.round(p.electrode_center_region_area_um2))} μm²` : '待计算';
  $('legalCorridors').textContent = isPad ? `${fmt(r?.candidate_library?.complete_pad_columns_after_pruning)} 条` : isVector ? `${fmt(r?.candidate_library?.path_columns)} 条已验证候选路线` : p ? `${fmt(p.vector_legal_single_track_corridors)} / ${fmt(g.graph_edges)}` : '待检查';
  $('junctionWindows').textContent = p ? fmt(p.ordered_junction_windows) : '待计算';
  $('topologyMatch').textContent = !activeJob && g.full_graph_matches_vector_topology ? '通过' : '待验证';
  $('scopeText').textContent = activeJob && isPad ?
    '结果边界：正在构造完整电极—Pad 网络；路径提案或内部认证数量尚不是最终电极数。导出 GDS 回读与独立整数审计通过后，才报告可行下界 L。' : isPad ?
    `结果边界：L=${fmt(cap?.lower_bound)} 是本次导出 GDS 经回读和整数几何审计确认的电极—专属 Pad 完整网络数；若严格审计未通过，L 暂不赋值。容量上界 ${fmt(cap?.upper_bound)} 来自原结构的整数中心装填界与已验证外侧 Pad 的几何割。${cap?.declared_geometric_model_optimality_proven ? 'L=U，当前输入及规则下的几何模型已证明最多；完整制造 DRC 不在证明范围。' : 'L 与 U 尚未闭合；有限候选库的最优值不能证明连续问题最多。'}` :
    '结果边界：几何试布线模式只认证电极到几何终端，不能当作 Pad 连通数量；完整 Pad 模式需另行运行。';
}

async function selectFile(id, spec = null) {
  state.detailAbort?.abort(); state.detailAbort = new AbortController();
  const file = state.files.find(f => f.id === id); if (!file) return;
  const task = taskForFile(file);
  if (!spec && task?.support_layer) spec = task.support_layer;
  state.selected = file; state.detail = null; state.view = 'support';
  state.job = task || null;
  state.jobFullId = null;
  state.jobFullStatus = null;
  localStorage.setItem('gds-workbench-file', id);
  renderFiles(); setStateBadge('读取中…');
  $('selectedTitle').textContent = file.name;
  $('selectedSubtitle').textContent = file.relative_path;
  $('detailStatus').textContent = '读取中';
  ['holeMetric', 'componentMetric', 'edgeMetric', 'outletMetric'].forEach(k => $(k).textContent = '—');
  ['pathLabel', 'sizeLabel', 'hashLabel', 'layerCount'].forEach(k => $(k).textContent = '—');
  $('layerList').replaceChildren(create('p', 'muted-copy', '正在读取 GDS 图层…'));
  $('layerReason').textContent = '读取中…';
  $('supportLayer').replaceChildren();
  $('copyHash').disabled = true;
  $('runButton').disabled = true; $('supportLayer').disabled = true;
  renderView(); renderAnalysis(); renderJob();
  const token = ++state.detailsToken;
  try {
    const [detail, fullJob] = await Promise.all([
      getJSON(`/api/file/${id}${spec ? `?spec=${encodeURIComponent(spec)}` : ''}`,
        { signal: AbortSignal.any([state.detailAbort.signal, AbortSignal.timeout(15000)]) }),
      task ? getJSON(`/api/job/${task.id}?detail=preview`,
        { cache: 'default', signal: AbortSignal.any([state.detailAbort.signal, AbortSignal.timeout(15000)]) }).catch(() => null) : Promise.resolve(null)
    ]);
    if (token !== state.detailsToken) return;
    state.detail = detail;
    if (fullJob) {
      state.job = fullJob;
      state.jobFullId = fullJob.id;
      state.jobFullStatus = fullJob.status;
      const previous = fullJob.rules || fullJob.result?.process_geometry?.rules || {};
      for (const [key, id] of [['electrode_diameter_um','electrodeDiameter'],['wire_width_um','wireWidth'],['spacing_um','spacing'],['margin_um','margin'],['collector_clearance_um','collectorClearance']]) {
        if (previous[key] != null) $(id).value = String(previous[key]);
      }
      $('minimumCenterSpacingMm').value = String((previous.minimum_center_spacing_um ?? 70) / 1000);
      $('electrodeRegionRadiusMm').value = String((previous.electrode_region_radius_um ?? 3000) / 1000);
      const previousPad = fullJob.pad_settings || fullJob.result?.routing?.pad_settings || {};
      for (const [key, id, fallback] of [['square_side_um','padFrameSideMm',32000],['pad_width_um','padWidthMm',500],['pad_length_um','padLengthMm',3000],['pad_pitch_um','padPitchMm',1000]])
        $(id).value = String((previousPad[key] ?? fallback) / 1000);
      $('padSizingMode').value = previousPad.sizing_mode || 'expand_frame';
      $('minimumPadWidthMm').value = String((previousPad.minimum_pad_width_um ?? previousPad.pad_width_um ?? 500) / 1000);
      $('minimumPadLengthMm').value = String((previousPad.minimum_pad_length_um ?? previousPad.pad_length_um ?? 3000) / 1000);
      if (fullJob.outlet_mode) $('outletMode').value = fullJob.outlet_mode;
      if (fullJob.method && ['four_side_pads','attached_islands'].includes(fullJob.method)) $('routingMode').value = fullJob.method;
    } else {
      $('collectorClearance').value = '25';
      $('minimumCenterSpacingMm').value = '0.07';
      $('electrodeRegionRadiusMm').value = '3';
      $('outletMode').value = 'auto_geometry';
      $('routingMode').value = 'four_side_pads';
      for (const [id, value] of [['padFrameSideMm',32],['padWidthMm',0.5],['padLengthMm',3],['padPitchMm',1]]) $(id).value = String(value);
      $('padSizingMode').value = 'expand_frame';
      $('minimumPadWidthMm').value = '0.5'; $('minimumPadLengthMm').value = '3';
    }
    if (!fullJob?.outlet_mode) $('outletMode').value = 'auto_geometry';
    renderDetails(); renderView(); renderJob();
  } catch (error) { if (token === state.detailsToken) { setStateBadge('读取失败', 'warn'); toast(`GDS 读取失败：${error.message}`, true); } }
}

function rulesPayload() {
  const fields = [['electrode_diameter_um', 'electrodeDiameter'], ['wire_width_um', 'wireWidth'], ['spacing_um', 'spacing'], ['margin_um', 'margin']];
  if ($('routingMode').value !== 'four_side_pads') fields.push(['collector_clearance_um', 'collectorClearance']);
  const rules = {};
  fields.forEach(([key, id]) => {
    const input = $(id), value = Number(input.value);
    const allowZero = key === 'spacing_um' || key === 'margin_um';
    const label = input.closest('label')?.firstChild.textContent?.trim() || id;
    if (input.value.trim() === '' || !Number.isFinite(value) || value < 0 || (!allowZero && value === 0))
      throw new Error(`${label} 必须是${allowZero ? '大于等于 0 的有限数值' : '大于 0 的有限数值'}`);
    rules[key] = value;
  });
  const centerInput = $('minimumCenterSpacingMm'), centerMm = Number(centerInput.value);
  if (centerInput.value.trim() === '' || !Number.isFinite(centerMm * 1000) || centerMm < 0)
    throw new Error('最小电极中心距必须是大于等于 0 的有限数值');
  rules.minimum_center_spacing_um = centerMm * 1000;
  const radiusInput = $('electrodeRegionRadiusMm'), radiusMm = Number(radiusInput.value);
  if (radiusInput.value.trim() === '' || !Number.isFinite(radiusMm * 1000) || radiusMm <= 0)
    throw new Error('电极中心允许半径必须是大于 0 的有限数值，单位 mm');
  rules.electrode_region_radius_um = radiusMm * 1000;
  return rules;
}

function padSettingsPayload(rules) {
  const fields = [['square_side_um','padFrameSideMm'],['pad_width_um','padWidthMm'],['pad_length_um','padLengthMm'],['pad_pitch_um','padPitchMm']];
  const values = {};
  for (const [key, id] of fields) {
    const input = $(id);
    const mm = Number(input.value);
    if (input.value.trim() === '' || !Number.isFinite(mm) || mm <= 0) throw new Error(`${id} 必须是正数`);
    values[key] = mm * 1000;
  }
  const contact = rules.wire_width_um + 2 * rules.margin_um;
  if (values.pad_width_um < contact || values.pad_length_um < contact)
    throw new Error(`Pad 宽度和长度至少需 ${contact} μm，才能容纳导线及两侧余量`);
  if (values.pad_pitch_um < values.pad_width_um + rules.spacing_um)
    throw new Error('Pad 节距必须不小于 Pad 宽度加网络间距');
  values.sizing_mode = $('padSizingMode').value;
  for (const [key, id, preferred] of [['minimum_pad_width_um','minimumPadWidthMm','pad_width_um'],['minimum_pad_length_um','minimumPadLengthMm','pad_length_um']]) {
    const dimension = Number($(id).value) * 1000;
    if ($(id).value.trim() === '' || !Number.isFinite(dimension) || dimension < contact || dimension > values[preferred])
      throw new Error(`Pad 最小尺寸须在 ${contact} μm 与参考尺寸之间`);
    values[key] = dimension;
  }
  return values;
}

async function startJob() {
  if (!state.detail || !state.detail.routable_input) return;
  if (state.submitPending || ['queued', 'running'].includes(taskForFile(state.detail)?.status)) return;
  let rules, padSettings;
  try {
    rules = rulesPayload();
    if ($('routingMode').value === 'four_side_pads') padSettings = padSettingsPayload(rules);
  } catch (error) { return toast(error.message, true); }
  state.submitPending = true;
  $('runButton').disabled = true;
  try {
    const data = await getJSON('/api/jobs', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ file_id: state.detail.id, support_layer: state.detail.selected_support_layer, rules, pad_settings: padSettings, outlet_mode: $('outletMode').value, method: $('routingMode').value }) });
    state.job = { id: data.job_id, file_id: state.detail.id, support_layer: state.detail.selected_support_layer, method: $('routingMode').value, pad_settings: padSettings, status: 'queued', stage: '排队中', progress: 0, message: '正在准备分析', artifacts: [] };
    state.jobFullId = null;
    state.jobFullStatus = null;
    state.jobsByFile[state.detail.id] = state.job;
    state.view = 'support';
    renderFiles(); renderJob(); renderDetails(); renderView();
    pollQueue();
  } catch (error) { toast(`无法启动：${error.message}`, true); }
  finally { state.submitPending = false; renderDetails(); }
}

async function pollQueue() {
  if (state.pollBusy || document.hidden) return;
  state.pollBusy = true;
  try {
    const previousSignature = fileStatusSignature();
    const previousSelected = `${state.job?.id || ''}:${state.job?.status || ''}`;
    const previousRender = jobRenderKey();
    const queue = await getJSON('/api/queue');
    $('queueIndicator').textContent = `运行 ${queue.running}/${queue.max_workers} · 排队 ${queue.queued}`;
    for (const [fileId, task] of Object.entries(queue.jobs || {})) {
      const old = state.jobsByFile[fileId];
      if (!old || !old.created_at || task.created_at >= old.created_at) state.jobsByFile[fileId] = task;
    }
    const task = taskForFile(state.selected);
    if (task && state.detail && task.id !== state.job?.id) {
      if (task.support_layer && task.support_layer !== state.detail.selected_support_layer) {
        await selectFile(state.selected.id, task.support_layer);
      } else { state.job = task; state.jobFullId = null; state.jobFullStatus = null; }
    } else if (task && state.job?.id === task.id) {
      state.job = { ...state.job, ...task };
    }
    if (task && state.job?.id === task.id &&
        ['complete', 'error', 'interrupted'].includes(task.status) &&
        (state.jobFullId !== task.id || state.jobFullStatus !== task.status)) {
      const currentId = task.id;
      const full = await getJSON(`/api/job/${currentId}?detail=preview`, { cache: 'default' });
      if (state.selected?.id === full.file_id && state.job?.id === currentId) {
        state.job = full; state.jobFullId = currentId; state.jobFullStatus = full.status;
        if (task.status === 'complete') {
          toast(full.message || '分析完成');
          state.view = full.artifacts?.includes('routing.gds') ? 'routing' : 'graph';
          renderView();
        } else if (task.status === 'error') toast(full.message || '运行失败', true);
      }
    }
    state.queueErrorShown = false;
    if (fileStatusSignature() !== previousSignature) renderFiles();
    if (jobRenderKey() !== previousRender) {
      renderJob();
      if (`${state.job?.id || ''}:${state.job?.status || ''}` !== previousSelected) renderDetails();
      else renderAnalysis();
      renderStages();
    }
  } catch (error) { if (!state.queueErrorShown) toast(`无法获取任务状态：${error.message}`, true); state.queueErrorShown = true; }
  finally { state.pollBusy = false; }
}

function jobRenderKey() {
  const job = state.job;
  return [job?.id, job?.status, job?.progress, job?.stage, job?.message,
    (job?.artifacts || []).join(','), state.jobFullId, state.jobFullStatus, historicalResultNote()].join('|');
}

function renderJob() {
  const job = jobMatches() ? state.job : null;
  updateElectrodeRegion();
  renderOverview(job);
  window.DistributionPanel?.update(job, state.detail);
  $('resultEmpty').hidden = !!job; $('resultContent').hidden = !job;
  if (!job) { $('resultBadge').textContent = '尚未运行'; $('resultMetrics').hidden = true; $('resultMetrics').style.display = 'none'; $('curveInfo').hidden = true; return; }
  $('resultBadge').textContent = job.status === 'complete' && padRun(job.result) && job.result.routing?.pad_connection_verified ? '电极—Pad 已回读验证' :
    job.status === 'complete' && padRun(job.result) ? '未构造完整 Pad 网络' :
    job.status === 'complete' && !job.artifacts?.includes('routing.gds') ? '仅几何分析' :
    job.status === 'complete' && job.result?.method === 'attached_islands' ? '几何试布线 · 未接 Pad' :
    ({ queued: '排队中', running: '运行中', complete: '已完成', error: '失败', interrupted: '已中断' })[job.status] || job.status;
  $('jobStage').textContent = job.stage || '准备中';
  $('jobPercent').textContent = `${Math.round((job.progress || 0) * 100)}%`;
  $('progressFill').style.width = `${Math.round((job.progress || 0) * 100)}%`;
  const excluded = job.result?.geometry?.navigation_geometry_certificate?.excluded_component_count;
  $('jobMessage').textContent = `${job.message || '—'}${job.status === 'complete' && job.result?.method === 'attached_islands' ? `；${terminalSummary(job.result)}；尚未验证到实体 Pad 的连接` : ''}${excluded ? `；导航域排除 ${excluded} 个不能满足本次几何终端间隔的分量，详见诊断；Pad 模式仍逐条验证完整连通` : ''}`;
  if (historicalResultNote()) $('jobMessage').textContent += `；${historicalResultNote()}`;
  const applied = job.rules || job.result?.process_geometry?.rules;
  const appliedPad = job.pad_settings || job.result?.routing?.pad_settings;
  $('jobRules').textContent = applied ? `本次规则：电极 ${applied.electrode_diameter_um} μm · 中心距 ${((applied.minimum_center_spacing_um || 0) / 1000).toFixed(3)} mm · 线宽 ${applied.wire_width_um} · 网络间距 ${applied.spacing_um} · 边缘余量 ${applied.margin_um} μm${appliedPad && job.method === 'four_side_pads' ? ` · Pad ${(appliedPad.pad_width_um / 1000).toFixed(3)} × ${(appliedPad.pad_length_um / 1000).toFixed(3)} mm · 节距 ${(appliedPad.pad_pitch_um / 1000).toFixed(3)} mm · 起始外框 ${(appliedPad.square_side_um / 1000).toFixed(3)} mm` : ''}` : '—';
  const result = job.result;
  if (applied) $('jobRules').textContent += applied.electrode_region_radius_um == null ?
    ' · 电极区域半径未限定（历史任务）' : ` · 电极区域半径 ${fmt(applied.electrode_region_radius_um / 1000)} mm`;
  $('resultMetrics').hidden = !result;
  $('resultMetrics').style.display = result ? 'grid' : 'none';
  if (result) {
    $('lowerBoundLabel').textContent = padRun(result) ? '已接通 Pad 下界 L' : '几何终端试布线下界 L';
    $('upperBoundLabel').textContent = padRun(result) ? '连续几何上界 U' : '当前锚点域几何上界 U';
    $('lowerBound').textContent = result.capacity_interval?.lower_bound == null ? '未定义' : String(result.capacity_interval.lower_bound);
    $('upperBound').textContent = result.capacity_interval?.upper_bound == null ? '未证明' : String(result.capacity_interval.upper_bound);
    $('roundtripStatus').textContent = result.routing?.gds_roundtrip_audit?.passed ? (padRun(result) ? 'Pad 连通通过' : '几何端通过') : '未运行';
  }
  const metrics = result?.routing?.preview_metrics;
  const curves = result?.routing?.routes?.map(route => route.curve).filter(Boolean) || [];
  const curveCount = metrics?.curve_count ?? curves.length;
  const corners = metrics?.curved_corner_count ?? curves.reduce((sum, curve) => sum + (curve.curved_corner_count || 0), 0);
  const radii = curves.map(curve => curve.minimum_bend_radius_um).filter(value => value != null);
  const minimumRadius = metrics?.minimum_bend_radius_um ?? (radii.length ? Math.min(...radii) : null);
  $('curveInfo').hidden = !curveCount;
  if (curveCount) $('curveInfo').textContent = `解析曲线：${fmt(corners)} 个转角采用 G¹ 切向连续曲线；最小数学弯曲半径 ${minimumRadius == null ? '—' : minimumRadius.toFixed(2)} μm。金属由弦偏差 ≤0.01 μm 的采样中心线生成，GDS 多边形另经坐标量化和回读核查。此半径尚未作为工艺下限认证。`;
  const links = $('artifactList'); links.replaceChildren();
  const names = { 'summary.json': '完整报告 JSON', 'navigation_diagnostics.json': '导航余量与孤岛证明 JSON', 'graph.png': '结构提取图 PNG', 'routing.png': '试布线预览 PNG', 'routing.gds': '金属 GDS', 'graph.json.gz': '走廊图数据', 'regions.json.gz': '矢量可放区' };
  (job.artifacts || []).forEach(file => {
    const a = create('a', '', names[file] || file);
    a.href = `/api/artifact/${job.id}/${encodeURIComponent(file)}`;
    if (!file.endsWith('.png')) a.download = file;
    links.append(a);
  });
}

function renderOverview(job) {
  const result = job?.result;
  const auditPassed = result?.routing?.gds_roundtrip_audit?.passed === true;
  const routes = result?.routing?.retained_routes;
  const count = auditPassed ? routes : result && routes === 0 ? 0 : null;
  $('electrodeHeroCount').textContent = count == null ? (job && ['queued', 'running'].includes(job.status) ? '…' : '—') : fmt(count);
  $('electrodeHeroDetail').textContent = auditPassed ?
    (padRun(result) ? `个电极分别连到专属 Pad，完整金属及承载桥已通过 GDS 回读。${result.capacity_interval?.integer_polygon_lower_verified ? '严格整数几何审计也已通过。' : '严格整数几何审计未通过，暂不计入证明下界。'}${result.capacity_interval?.declared_geometric_model_optimality_proven ? 'L=U，已证明当前几何模型内最多。' : '当前尚未证明最多。'}本次外框可用槽位 ${fmt(result.routing.pad_slot_count)} 个，并非固定总上限；只导出接通的 Pad。` :
    `个电极接到几何候选点并通过 GDS 回读。${terminalSummary(result)}；尚未接到实体 Pad。`) :
    result && routes === 0 ? '本次参数下没有形成完整合法网络；0 是本次构造结果，不是不可布线证明。' :
    job && ['queued', 'running'].includes(job.status) ? '正在运行；完成实体金属导出和回读后显示布置数量。' :
    '选择 case 并运行后显示数量；左侧已完成的 case 可直接查看历史结果。';
  const rules = job?.rules || result?.routing?.rules;
  const spacing = rules?.minimum_center_spacing_um;
  $('overviewSpacing').textContent = spacing == null ? '—' : spacing === 0 ? '未设置' : `${fmt(spacing / 1000)} mm`;
  const exitPolicy = result?.routing?.outer_exit_policy;
  $('overviewExitWindow').textContent = exitPolicy?.maximum_launch_depth_from_envelope_um != null ? `外包络内 ${exitPolicy.maximum_launch_depth_from_envelope_um.toFixed(3)} μm` : exitPolicy ? `${(exitPolicy.minimum_launch_radius_um / 1000).toFixed(4)}–${(exitPolicy.source_max_radius_um / 1000).toFixed(4)} mm` : '待提取';
  $('overviewUpper').textContent = result?.capacity_interval?.upper_bound == null ? '—' : fmt(result.capacity_interval.upper_bound);
  $('overviewAudit').textContent = auditPassed ? (result?.capacity_interval?.declared_geometric_model_optimality_proven ? '模型内最多已证明' : padRun(result) ? result.routing.integer_polygon_audit?.outer_exit_policy_verified ? 'Pad 与外端出口通过' : 'Pad 连通通过' : '几何端通过') : result && routes === 0 ? '无 GDS 输出' : job?.status === 'error' ? '运行失败' : '待完成';
  $('overviewAudit').classList.toggle('passed', auditPassed);
}

$('zoomIn').addEventListener('click', () => zoom(1.8));
function updateElectrodeRegion() {
  const value = Number($('electrodeRegionRadiusMm').value);
  const reference = state.detail?.geometry?.electrode_region_reference_um;
  const draft = Number.isFinite(value) && value > 0 && reference ?
    { radius_um: value * 1000, reference_um: reference } : null;
  const job = jobMatches() ? state.job : null;
  const actual = job?.result?.routing?.electrode_region;
  const shown = state.view === 'support' ? draft : actual?.enabled ? actual : null;
  if (state.viewer.kind === 'vector') vectorPreview?.setRegion(shown);
  if (state.viewer.kind === 'svg') {
    state.viewer.svg.querySelector('[data-electrode-region]')?.remove();
    if (shown) {
      const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
      for (const [key, value] of Object.entries({cx:shown.reference_um[0],cy:shown.reference_um[1],r:shown.radius_um,
        fill:'none',stroke:'#c56e19','stroke-width':'1.5','stroke-dasharray':'7 5','vector-effect':'non-scaling-stroke','data-electrode-region':'true'}))
        circle.setAttribute(key,value);
      state.viewer.svg.querySelector('g')?.append(circle);
    }
  }
  const applied = job?.rules?.electrode_region_radius_um;
  $('electrodeRegionResult').textContent = actual?.enabled ?
    `当前结果使用半径 ${fmt(actual.radius_um / 1000)} mm；最远电极 ${actual.maximum_selected_radius_um == null ? '无电极' : `${fmt(actual.maximum_selected_radius_um / 1000)} mm`}。原始支撑页橙色虚线预览下次运行范围。` :
    job && applied == null ? '当前显示的是未限制半径的历史任务；下次运行使用上方新参数。原始支撑页橙色虚线可预览范围。' :
    `下次运行：半径 ${Number.isFinite(value) && value > 0 ? fmt(value) : '—'} mm，直径 ${Number.isFinite(value) && value > 0 ? fmt(2 * value) : '—'} mm。原始支撑页橙色虚线为允许范围。`;
}
$('electrodeRegionRadiusMm').addEventListener('input', updateElectrodeRegion);
$('zoomOut').addEventListener('click', () => zoom(1 / 1.8));
$('zoomReset').addEventListener('click', resetZoom);
const canvasStage = $('canvasStage');
canvasStage.addEventListener('wheel', event => {
  if (!state.viewer.kind) return;
  event.preventDefault();
  zoom(event.deltaY < 0 ? 1.25 : 1 / 1.25, event.clientX, event.clientY);
}, { passive: false });
canvasStage.addEventListener('pointerdown', event => {
  if (!state.viewer.kind || !event.isPrimary || event.button !== 0 ||
      event.target.closest('.zoom-controls') || state.viewer.drag) return;
  event.preventDefault();
  canvasStage.setPointerCapture(event.pointerId);
  state.viewer.drag = { pointerId: event.pointerId, x: event.clientX, y: event.clientY };
  canvasStage.classList.add('dragging');
});
canvasStage.addEventListener('pointermove', event => {
  const v = state.viewer;
  if (v.drag?.pointerId !== event.pointerId) return;
  if (event.pointerType === 'mouse' && event.buttons === 0) { cancelDrag(event.pointerId); return; }
  if (v.kind === 'vector') {
    vectorPreview.pan(event.clientX - v.drag.x, event.clientY - v.drag.y);
  } else if (v.kind === 'svg') {
    const before = svgPoint(v.drag.x, v.drag.y), now = svgPoint(event.clientX, event.clientY);
    const box = v.svg.viewBox.baseVal;
    v.svg.setAttribute('viewBox', `${box.x + before.x - now.x} ${box.y + before.y - now.y} ${box.width} ${box.height}`);
  } else if (v.kind === 'image') {
    v.imageX += event.clientX - v.drag.x; v.imageY += event.clientY - v.drag.y;
    $('mainImage').style.transform = `translate(${v.imageX}px, ${v.imageY}px) scale(${v.imageScale})`;
  }
  v.drag.x = event.clientX; v.drag.y = event.clientY;
});
canvasStage.addEventListener('pointerup', event => cancelDrag(event.pointerId));
canvasStage.addEventListener('pointercancel', event => cancelDrag(event.pointerId));
canvasStage.addEventListener('lostpointercapture', event => cancelDrag(event.pointerId));
canvasStage.addEventListener('dragstart', event => {
  if (!event.target.closest('.zoom-controls')) event.preventDefault();
});
$('refreshButton').addEventListener('click', refreshFiles);
$('searchInput').addEventListener('input', renderFiles);
$('runButton').addEventListener('click', startJob);
$('outletMode').addEventListener('change', renderDetails);
$('routingMode').addEventListener('change', renderDetails);
$('supportLayer').addEventListener('change', () => { if (state.selected) selectFile(state.selected.id, $('supportLayer').value); });
$('copyHash').addEventListener('click', async () => { try { await navigator.clipboard.writeText(state.detail?.sha256 || ''); toast('SHA-256 已复制'); } catch { toast('复制失败', true); } });
document.addEventListener('keydown', event => { if (event.key === '/' && !['INPUT', 'TEXTAREA'].includes(document.activeElement.tagName)) { event.preventDefault(); $('searchInput').focus(); } });
document.addEventListener('visibilitychange', () => { if (!document.hidden) pollQueue(); });
refreshFiles();
