/* Exact vector commands rendered off the UI thread, with viewport culling. */
const scenes = new Map(), loading = new Map(), fetches = new Map(), thumbnails = [];
let sceneBytes = 0, mainURL = null, pendingMain = null, pumping = false, activeLoad = 0;
const MAX_BYTES = 160 * 1024 * 1024;
const pause = () => new Promise(resolve => setTimeout(resolve, 0));

async function sceneFor(url) {
  if (scenes.has(url)) {
    const scene = scenes.get(url); scenes.delete(url); scenes.set(url, scene); return scene;
  }
  if (loading.has(url)) return loading.get(url);
  const controller = new AbortController(); fetches.set(url, controller);
  const promise = (async () => {
    const response = await fetch(url, { signal: controller.signal });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const data = await response.json(), started = performance.now();
    const paths = [], bins = new Map(), large = [], base = data.viewBox;
    const world = [base[0], data.flipY - base[1] - base[3], base[2], base[3]];
    const cellX = Math.max(world[2] / 48, 1e-9), cellY = Math.max(world[3] / 48, 1e-9);
    let bytes = 0;
    for (let index = 0; index < data.paths.length; index++) {
      const record = data.paths[index]; bytes += record.d.length * 2 + 128;
      paths.push({ ...record, d: undefined, path: new Path2D(record.d), index });
      const b = record.bounds;
      const x1 = Math.floor((b[0] - world[0]) / cellX), x2 = Math.floor((b[2] - world[0]) / cellX);
      const y1 = Math.floor((b[1] - world[1]) / cellY), y2 = Math.floor((b[3] - world[1]) / cellY);
      if ((x2 - x1 + 1) * (y2 - y1 + 1) > 64) large.push(index);
      else for (let x = x1; x <= x2; x++) for (let y = y1; y <= y2; y++) {
        const key = `${x}:${y}`; if (!bins.has(key)) bins.set(key, []); bins.get(key).push(index);
      }
      if (index % 160 === 159) {
        await pause(); controller.signal.throwIfAborted();
      }
    }
    controller.signal.throwIfAborted();
    const scene = { base, flipY: data.flipY, focus: data.focusBox, paths, bins, large,
      world, cellX, cellY, bytes, prepareMs: performance.now() - started };
    scenes.set(url, scene); sceneBytes += bytes;
    for (const [oldURL, old] of scenes) {
      if (sceneBytes <= MAX_BYTES) break;
      if (oldURL !== mainURL && oldURL !== url) { scenes.delete(oldURL); sceneBytes -= old.bytes; }
    }
    return scene;
  })();
  loading.set(url, promise);
  try { return await promise; } finally {
    if (loading.get(url) === promise) loading.delete(url);
    if (fetches.get(url) === controller) fetches.delete(url);
  }
}

function cancelLoading() {
  for (const [url, controller] of fetches) {
    controller.abort(); loading.delete(url);
  }
  fetches.clear();
}

function visiblePaths(scene, box, strokePad = 0) {
  const bounds = [box[0] - strokePad, scene.flipY - box[1] - box[3] - strokePad,
    box[0] + box[2] + strokePad, scene.flipY - box[1] + strokePad];
  if (bounds[0] <= scene.world[0] && bounds[1] <= scene.world[1] &&
      bounds[2] >= scene.world[0] + scene.world[2] && bounds[3] >= scene.world[1] + scene.world[3]) return scene.paths;
  const ids = new Set(scene.large);
  const x1 = Math.max(-1, Math.floor((bounds[0] - scene.world[0]) / scene.cellX));
  const x2 = Math.min(49, Math.floor((bounds[2] - scene.world[0]) / scene.cellX));
  const y1 = Math.max(-1, Math.floor((bounds[1] - scene.world[1]) / scene.cellY));
  const y2 = Math.min(49, Math.floor((bounds[3] - scene.world[1]) / scene.cellY));
  for (let x = x1; x <= x2; x++) for (let y = y1; y <= y2; y++)
    for (const id of scene.bins.get(`${x}:${y}`) || []) ids.add(id);
  return [...ids].sort((a, b) => a - b).map(id => scene.paths[id]).filter(p =>
    p.bounds[0] <= bounds[2] && p.bounds[2] >= bounds[0] && p.bounds[1] <= bounds[3] && p.bounds[3] >= bounds[1]);
}

async function draw(request, scene) {
  const started = performance.now(), box = request.box || (request.focus && scene.focus) || scene.base;
  const width = Math.max(1, Math.round(request.width * request.dpr)), height = Math.max(1, Math.round(request.height * request.dpr));
  const canvas = new OffscreenCanvas(width, height), ctx = canvas.getContext('2d', { alpha: false });
  ctx.fillStyle = '#f7f8f5'; ctx.fillRect(0, 0, width, height);
  const scale = Math.min(width / box[2], height / box[3]);
  const tx = (width - box[2] * scale) / 2 - box[0] * scale;
  const ty = (height - box[3] * scale) / 2 - box[1] * scale;
  ctx.setTransform(scale, 0, 0, -scale, tx, ty + scene.flipY * scale);
  // preserveAspectRatio="xMidYMid meet" exposes extra world area when the
  // canvas aspect differs from the viewBox. Cull against that actual area.
  const viewportBox = [box[0] - (width / scale - box[2]) / 2,
    box[1] - (height / scale - box[3]) / 2, width / scale, height / scale];
  const edgePad = 1 / scale + (request.highlight ? 27.5 : request.focus ? 5 : 0);
  let visible = visiblePaths(scene, viewportBox, edgePad);
  if (request.highlight) visible = [...visible].sort((a, b) => order(a) - order(b) || a.index - b.index);
  for (let index = 0; index < visible.length; index++) {
    const record = visible[index];
    ctx.globalAlpha = request.highlight && record.fill === '#3b807d' ? .24 : record.opacity;
    if (record.fill !== 'none') { ctx.fillStyle = record.fill; ctx.fill(record.path, 'evenodd'); }
    const highlight = request.highlight && ['#da683e', '#142e51'].includes(record.fill);
    if (record.stroke || highlight) {
      ctx.strokeStyle = record.stroke || record.fill;
      ctx.lineWidth = highlight ? 55 : request.focus && record.fill === 'none' ? 10 : record.width;
      ctx.lineJoin = 'round'; ctx.lineCap = 'round'; ctx.stroke(record.path);
    }
    if (index % 250 === 249) {
      await pause();
      if (request.role === 'main' && (pendingMain || request.token !== activeLoad)) return;
      if (request.role !== 'main' && pendingMain && request.token === activeLoad) {
        thumbnails.unshift(request); return;
      }
    }
  }
  if (request.token !== activeLoad) return;
  const bitmap = canvas.transferToImageBitmap();
  self.postMessage({ ...request, type: 'frame', box, bitmap, visible: visible.length,
    total: scene.paths.length, renderMs: performance.now() - started }, [bitmap]);
}

function order(record) {
  return record.fill === '#3b807d' ? 0 : record.fill === '#da683e' ? 1 : record.fill === 'none' ? 3 : 2;
}

async function pump() {
  if (pumping) return; pumping = true;
  try {
    while (pendingMain || thumbnails.length) {
      const request = pendingMain || thumbnails.shift();
      if (request.role === 'main') pendingMain = null;
      try {
        const scene = await sceneFor(request.url);
        if (request.token !== activeLoad) continue;
        await draw(request, scene);
      } catch (error) { self.postMessage({ type: 'error', id: request.id, token: request.token, message: error.message }); }
    }
  } finally { pumping = false; }
}

self.onmessage = event => {
  const request = event.data;
  if (request.type === 'load') {
    mainURL = request.url; activeLoad = request.token; pendingMain = null;
    // Obsolete case thumbnails should never delay the newly selected case.
    thumbnails.length = 0;
    sceneFor(request.url).then(scene => {
      if (request.token === activeLoad) self.postMessage({ type: 'loaded', token: request.token,
        base: scene.base, flipY: scene.flipY, total: scene.paths.length, prepareMs: scene.prepareMs });
    }).catch(error => self.postMessage({ type: 'error', token: request.token, message: error.message }));
  } else if (request.type === 'render') {
    if (request.role === 'main') pendingMain = request;
    else thumbnails.push(request);
    pump();
  } else if (request.type === 'cancel') {
    activeLoad = request.token; pendingMain = null; thumbnails.length = 0; cancelLoading();
  }
};
