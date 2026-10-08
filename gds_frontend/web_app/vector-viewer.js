/* Retain vector geometry in a worker; move its latest frame immediately via CSS. */
class VectorPreview {
  constructor(stage, viewport) {
    this.stage = stage; this.viewport = viewport; this.token = 0; this.frame = 0;
    this.worker = new Worker('/vector-worker.js'); this.thumbnails = new Map();
    this.worker.onmessage = event => this.receive(event.data);
    this.worker.onerror = event => this.fail(new Error(event.message || '后台矢量绘制失败'));
    this.observer = new ResizeObserver(() => { if (this.base) this.schedule(); }); this.observer.observe(stage);
  }
  load(url, label) {
    this.cancel(); this.url = url; this.base = null; this.box = null; this.presented = null;
    this.started = performance.now();
    this.canvas = document.createElement('canvas'); this.canvas.className = 'vector-canvas';
    this.canvas.setAttribute('role', 'img'); this.canvas.setAttribute('aria-label', label);
    this.canvas.dataset.state = 'loading'; this.viewport.replaceChildren(this.canvas);
    this.overlay = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    this.overlay.style.cssText = 'position:absolute;inset:0;pointer-events:none;overflow:hidden';
    this.overlay.setAttribute('aria-label', '电极中心允许区域');
    this.viewport.append(this.overlay);
    this.context = this.canvas.getContext('bitmaprenderer');
    if (!this.context) this.context = this.canvas.getContext('2d');
    const token = this.token;
    this.worker.postMessage({ type: 'load', url, token });
    return new Promise((resolve, reject) => { this.loading = { token, resolve, reject }; });
  }
  cancel() {
    ++this.token; if (this.frame) cancelAnimationFrame(this.frame); this.frame = 0;
    if (this.loading) { this.loading.reject(new DOMException('View replaced', 'AbortError')); this.loading = null; }
    this.worker.postMessage({ type: 'cancel', token: this.token });
    this.thumbnails.clear(); this.base = null;
  }
  fail(error) {
    if (this.canvas) this.canvas.dataset.state = 'failed';
    if (this.loading) { this.loading.reject(error); this.loading = null; }
  }
  receive(message) {
    if (message.type === 'loaded' && message.token === this.token) {
      this.canvas.dataset.prepareMs = message.prepareMs.toFixed(1);
      this.base = message.base; this.box = [...message.base]; this.flipY = message.flipY; this.schedule();
    } else if (message.type === 'frame') {
      if (message.role !== 'main') {
        const target = this.thumbnails.get(message.id);
        if (message.token === this.token && target?.isConnected) {
          const context = target.getContext('bitmaprenderer') || target.getContext('2d');
          target.width = message.bitmap.width; target.height = message.bitmap.height;
          if (context.transferFromImageBitmap) context.transferFromImageBitmap(message.bitmap);
          else { context.drawImage(message.bitmap, 0, 0); message.bitmap.close(); }
          target.closest('.stage-thumbnail')?.classList.add('loaded');
        } else message.bitmap.close();
        this.thumbnails.delete(message.id); return;
      }
      if (message.token !== this.token) { message.bitmap.close(); return; }
      this.canvas.width = message.bitmap.width; this.canvas.height = message.bitmap.height;
      this.canvas.style.width = `${message.width}px`; this.canvas.style.height = `${message.height}px`;
      if (this.context.transferFromImageBitmap) this.context.transferFromImageBitmap(message.bitmap);
      else { this.context.drawImage(message.bitmap, 0, 0); message.bitmap.close(); }
      this.presented = { box: message.box, width: message.width, height: message.height };
      this.canvas.dataset.state = 'ready'; this.canvas.dataset.visiblePaths = message.visible;
      this.canvas.dataset.totalPaths = message.total; this.canvas.dataset.renderMs = message.renderMs.toFixed(1);
      this.canvas.dataset.viewBox = this.box.join(' '); this.transform();
      if (this.loading) {
        this.canvas.dataset.firstFrameMs = (performance.now() - this.started).toFixed(1);
        this.loading.resolve(); this.loading = null;
      }
    } else if (message.type === 'error' && message.token === this.token) {
      if (message.id) {
        const target = this.thumbnails.get(message.id);
        const wrapper = target?.closest('.stage-thumbnail');
        wrapper?.classList.add('failed');
        const label = wrapper?.querySelector('.stage-thumb-loading');
        if (label) label.textContent = '预览暂不可用';
        this.thumbnails.delete(message.id);
      } else this.fail(new Error(message.message));
    }
  }
  fit(box, width, height) {
    const scale = Math.min(width / box[2], height / box[3]);
    return { scale, x: (width - box[2] * scale) / 2 - box[0] * scale,
      y: (height - box[3] * scale) / 2 - box[1] * scale };
  }
  transform() {
    this.drawRegion();
    if (!this.presented || !this.box) return;
    const before = this.fit(this.presented.box, this.presented.width, this.presented.height);
    const after = this.fit(this.box, this.stage.clientWidth, this.stage.clientHeight);
    const scale = after.scale / before.scale;
    const x = after.x - before.x * scale, y = after.y - before.y * scale;
    this.canvas.style.transform = `translate(${x}px,${y}px) scale(${scale})`;
    this.canvas.dataset.viewBox = this.box.join(' ');
  }
  setRegion(region) { this.region = region; this.drawRegion(); }
  drawRegion() {
    if (!this.overlay) return;
    this.overlay.replaceChildren();
    if (!this.box || !this.region || !Number.isFinite(this.flipY)) return;
    this.overlay.setAttribute('viewBox', this.box.join(' '));
    const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    circle.setAttribute('cx', this.region.reference_um[0]);
    circle.setAttribute('cy', this.flipY - this.region.reference_um[1]);
    circle.setAttribute('r', this.region.radius_um);
    circle.setAttribute('fill', 'none'); circle.setAttribute('stroke', '#c56e19');
    circle.setAttribute('stroke-width', '1.5'); circle.setAttribute('stroke-dasharray', '7 5');
    circle.setAttribute('vector-effect', 'non-scaling-stroke');
    const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
    title.textContent = `电极中心允许半径 ${this.region.radius_um / 1000} mm`;
    circle.append(title); this.overlay.append(circle);
  }
  schedule() {
    this.transform();
    if (this.frame || !this.box) return;
    this.frame = requestAnimationFrame(() => {
      this.frame = 0;
      if (!this.box) return;
      this.worker.postMessage({ type: 'render', role: 'main', url: this.url, token: this.token,
        box: [...this.box], width: this.stage.clientWidth, height: this.stage.clientHeight,
        dpr: Math.min(devicePixelRatio || 1, 2) });
    });
  }
  zoom(factor, clientX, clientY) {
    if (!this.box) return;
    const rect = this.stage.getBoundingClientRect();
    const width = this.stage.clientWidth, height = this.stage.clientHeight;
    const left = rect.left + this.stage.clientLeft, top = rect.top + this.stage.clientTop;
    const fit = this.fit(this.box, width, height);
    const px = ((clientX ?? left + width / 2) - left - fit.x) / fit.scale;
    const py = ((clientY ?? top + height / 2) - top - fit.y) / fit.scale;
    const worldWidth = Math.max(this.base[2] / 500, Math.min(this.base[2] * 2, this.box[2] / factor));
    const ratio = worldWidth / this.box[2];
    this.box = [px - (px - this.box[0]) * ratio, py - (py - this.box[1]) * ratio, worldWidth, this.box[3] * ratio];
    this.schedule();
  }
  pan(dx, dy) {
    if (!this.box) return;
    const fit = this.fit(this.box, this.stage.clientWidth, this.stage.clientHeight);
    this.box[0] -= dx / fit.scale; this.box[1] -= dy / fit.scale; this.schedule();
  }
  reset() { if (this.base) { this.box = [...this.base]; this.schedule(); } }
  thumbnail(url, canvas, options = {}) {
    const id = `thumb-${++VectorPreview.nextID}`;
    this.thumbnails.set(id, canvas);
    this.worker.postMessage({ type: 'render', role: 'thumbnail', id, url, token: this.token,
      width: 240, height: 140, dpr: 1, ...options });
  }
}
VectorPreview.nextID = 0;
