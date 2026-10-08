const $ = id => document.getElementById(id);
const admin = { files: [], tasks: {}, selected: new Set(), layers: {}, choices: {}, queue: null, submitting: false, polling: false };
let renderedJobsKey = null;

async function json(url, options = {}) {
  const response = await fetch(url, { cache: 'no-store', ...options });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `HTTP ${response.status}`);
  return data;
}

function node(tag, className, value) {
  const element = document.createElement(tag);
  if (className) element.className = className;
  if (value != null) element.textContent = value;
  return element;
}

function message(value, error = false) {
  $('submitMessage').textContent = value;
  $('submitMessage').classList.toggle('error', error);
}

function taskFor(file) {
  let task = admin.tasks[file.id];
  if (task?.stale_reason === 'pad_layout_algorithm_changed' &&
      ['center_preferred_pad_routing_v12_electrode_region',
       'center_preferred_pad_routing_v13_attachment_cells'].includes(task.solver_revision)) {
    task = { ...task, stale_reason: null,
      parameter_note: task.solver_revision === 'center_preferred_pad_routing_v13_attachment_cells' ?
        '历史布局采用限时搜索；新任务取消默认求解时间限制。原结果保留供查看，尚未按不限时策略重跑。' :
        '历史布局采用旧位置搜索；新任务使用圆岛附着自适应搜索。原结果保留供查看，尚未按新策略重跑。' };
  }
  return task?.input_sha256 && task.input_sha256 !== file.sha256 ? null : task;
}

function active(file) {
  return ['queued', 'running'].includes(taskFor(file)?.status);
}

function fileStatusSignature() {
  return admin.files.map(file => {
    const task = taskFor(file);
    return `${file.id}:${task?.id || ''}:${task?.status || ''}:${task?.lower_bound ?? ''}:${task?.stale_reason || ''}:${task?.parameter_note || ''}`;
  }).join('|');
}

function statusText(task) {
  if (!task) return '待运行';
  if (task.stale_reason) return '算法已更新 · 请重跑';
  if (task.status === 'complete') return (task.lower_bound == null ? '已完成 · 无认证数量' : `已完成 · ${task.lower_bound} 个电极`) + (task.parameter_note ? ' · 历史布局' : '');
  return ({ queued: '排队中', running: '运行中', error: '运行失败', interrupted: '已中断' })[task.status] || task.status;
}

function updateSelection() {
  $('selectedCount').textContent = `${admin.selected.size} 已选择`;
  $('submitBatch').disabled = !admin.selected.size || admin.submitting;
  $('submitBatch').textContent = admin.submitting ? '正在校验并入队…' : `批量提交 ${admin.selected.size} 个任务 ↗`;
}

function renderFiles() {
  const list = $('adminFileList'); list.replaceChildren();
  for (const file of admin.files) {
    if (active(file)) admin.selected.delete(file.id);
    const task = taskFor(file);
    const row = node('div', `admin-file-item${admin.selected.has(file.id) ? ' selected' : ''}${active(file) ? ' disabled' : ''}`);
    const checkbox = node('input'); checkbox.type = 'checkbox'; checkbox.checked = admin.selected.has(file.id);
    checkbox.disabled = active(file); checkbox.setAttribute('aria-label', `选择 ${file.name}`);
    checkbox.addEventListener('change', () => {
      if (checkbox.checked) admin.selected.add(file.id); else admin.selected.delete(file.id);
      row.classList.toggle('selected', checkbox.checked); updateSelection();
    });
    const copy = node('div', 'admin-file-copy');
    copy.append(node('strong', 'admin-file-name', file.name), node('span', 'admin-file-path', file.relative_path));
    copy.append(node('span', `admin-file-status ${task?.status || ''}`, statusText(task)));
    const spec = node('select'); spec.setAttribute('aria-label', `${file.name} 支撑层`);
    const autoOption = node('option', '', admin.layers[file.id] ? `自动建议 ${admin.layers[file.id].suggested_support_layer}` : '自动建议支撑层');
    autoOption.value = 'auto'; spec.append(autoOption);
    if (admin.layers[file.id]) {
      for (const [layer, count] of Object.entries(admin.layers[file.id].layers)) {
        const option = node('option', '', `${layer} · ${count} 多边形`);
        option.value = layer; spec.append(option);
      }
    }
    spec.value = admin.choices[file.id] || 'auto';
    spec.addEventListener('change', () => { admin.choices[file.id] = spec.value; });
    copy.append(spec);
    const load = node('button', 'admin-load-layer', admin.layers[file.id] ? '已读取图层' : '读取图层后手动指定');
    load.type = 'button'; load.disabled = !!admin.layers[file.id];
    load.addEventListener('click', async () => {
      load.disabled = true; load.textContent = '读取中…';
      try {
        admin.layers[file.id] = await json(`/api/file/${file.id}`);
        renderFiles();
      } catch (error) { load.disabled = false; load.textContent = '读取失败，重试'; message(`${file.name}: ${error.message}`, true); }
    });
    copy.append(load);
    row.append(checkbox, copy); list.append(row);
  }
  updateSelection();
}

function renderJobs() {
  const queue = admin.queue;
  const key = [queue?.running, queue?.queued, queue?.max_workers,
    ...admin.files.map(file => {
      const task = taskFor(file);
      return [file.id, task?.id, task?.status, task?.stage, task?.progress, task?.lower_bound,
        task?.stale_reason, (task?.artifacts || []).join(',')].join(':');
    })].join('|');
  if (key === renderedJobsKey) return;
  renderedJobsKey = key;
  $('queueTop').textContent = queue ? `运行 ${queue.running}/${queue.max_workers} · 排队 ${queue.queued}` : '读取队列…';
  $('queueSummary').textContent = queue ? `${queue.running} 运行中 · ${queue.queued} 排队 · 同时最多 ${queue.max_workers}` : '—';
  const list = $('adminJobs'); list.replaceChildren();
  const jobs = admin.files.map(file => ({ file, job: taskFor(file) })).filter(x => x.job);
  jobs.sort((a, b) => (b.job.created_at || '').localeCompare(a.job.created_at || ''));
  if (!jobs.length) { list.append(node('p', 'admin-empty', '目前没有任务记录。')); return; }
  for (const { file, job } of jobs) {
    const row = node('div', 'admin-job');
    const title = node('div'); title.append(node('strong', '', file.name));
    title.append(node('small', '', job.batch_id ? `批次 ${job.batch_id.slice(0, 8)} · ${job.support_layer || '自动层'}` : `单个任务 · ${job.support_layer || '自动层'}`));
    const status = node('div', `admin-job-status ${job.status || ''}`, statusText(job));
    const progress = node('div', 'admin-job-progress');
    progress.append(node('span', '', `${job.stage || statusText(job)} · ${Math.round((job.progress || 0) * 100)}%`));
    const track = node('div', 'progress-track');
    const fill = node('span'); fill.style.width = `${Math.round((job.progress || 0) * 100)}%`;
    track.append(fill); progress.append(track);
    const links = node('div', 'admin-job-links');
    const open = node('a', '', '打开 case ↗'); open.href = `/?file=${encodeURIComponent(file.id)}`; links.append(open);
    if (job.artifacts?.includes('routing.gds')) {
      const gds = node('a', '', 'GDS ↓'); gds.href = `/api/artifact/${job.id}/routing.gds`; links.append(gds);
    }
    if (job.artifacts?.includes('summary.json')) {
      const summary = node('a', '', '报告 ↓'); summary.href = `/api/artifact/${job.id}/summary.json`; links.append(summary);
    }
    row.append(title, status, progress, links); list.append(row);
  }
}

async function refresh() {
  try {
    const [files, queue] = await Promise.all([json('/api/files'), json('/api/queue')]);
    admin.files = files.files.filter(file => file.routable_input);
    admin.tasks = { ...(files.tasks || {}), ...(queue.jobs || {}) };
    admin.queue = queue;
    $('inputRoot').textContent = files.root;
    renderFiles(); renderJobs();
  } catch (error) { message(`刷新失败：${error.message}`, true); }
}

async function poll() {
  if (admin.polling || document.hidden) return;
  admin.polling = true;
  try {
    const previousSignature = fileStatusSignature();
    const queue = await json('/api/queue');
    admin.queue = queue;
    admin.tasks = { ...admin.tasks, ...(queue.jobs || {}) };
    if (fileStatusSignature() !== previousSignature) renderFiles();
    renderJobs();
  } catch (error) { message(`队列状态读取失败：${error.message}`, true); }
  finally { admin.polling = false; }
}

function number(id, allowZero = false) {
  const input = $(id), value = Number(input.value);
  if (input.value.trim() === '' || !Number.isFinite(value) || value < 0 || (!allowZero && value === 0))
    throw new Error(`${input.closest('label')?.firstChild.textContent?.trim() || id} 必须是${allowZero ? '大于等于 0 的有限数值' : '大于 0 的有限数值'}`);
  return value;
}

function payload() {
  const method = $('adminMode').value;
  const rules = {
    electrode_diameter_um: number('adminElectrode'),
    wire_width_um: number('adminWire'),
    spacing_um: number('adminSpacing', true),
    margin_um: number('adminMargin', true),
    minimum_center_spacing_um: number('adminCenterSpacing', true) * 1000,
    electrode_region_radius_um: number('adminRegionRadius') * 1000
  };
  const result = {
    file_ids: [...admin.selected], support_layers: {}, method, rules,
    outlet_mode: method === 'four_side_pads' ? 'external_boundary' : $('adminOutlet').value
  };
  for (const id of admin.selected) if (admin.choices[id] && admin.choices[id] !== 'auto') result.support_layers[id] = admin.choices[id];
  if (method === 'four_side_pads') {
    result.pad_settings = {
      square_side_um: number('adminFrame') * 1000,
      pad_width_um: number('adminPadWidth') * 1000,
      pad_length_um: number('adminPadLength') * 1000,
      pad_pitch_um: number('adminPadPitch') * 1000,
      sizing_mode: $('adminPadSizing').value,
      minimum_pad_width_um: number('adminMinPadWidth') * 1000,
      minimum_pad_length_um: number('adminMinPadLength') * 1000
    };
    const contact = rules.wire_width_um + 2 * rules.margin_um;
    if (result.pad_settings.pad_width_um < contact || result.pad_settings.pad_length_um < contact)
      throw new Error(`Pad 宽度与长度至少 ${contact} μm`);
    if (result.pad_settings.pad_pitch_um < result.pad_settings.pad_width_um + rules.spacing_um)
      throw new Error('Pad 节距必须不小于 Pad 宽度加网络间距');
    if (result.pad_settings.minimum_pad_width_um < contact || result.pad_settings.minimum_pad_width_um > result.pad_settings.pad_width_um ||
        result.pad_settings.minimum_pad_length_um < contact || result.pad_settings.minimum_pad_length_um > result.pad_settings.pad_length_um)
      throw new Error(`Pad 最小尺寸须在 ${contact} μm 与参考尺寸之间`);
  } else rules.collector_clearance_um = number('adminCollector');
  return result;
}

async function submit() {
  if (!admin.selected.size || admin.submitting) return;
  let data;
  try { data = payload(); } catch (error) { message(error.message, true); return; }
  admin.submitting = true; updateSelection();
  try {
    const response = await json('/api/jobs/batch', {
      method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(data)
    });
    admin.selected.clear();
    message(`已提交 ${response.accepted} 个独立任务，批次 ${response.batch_id.slice(0, 8)}。队列将自动开始运行。`);
    await refresh();
  } catch (error) { message(`批量提交失败：${error.message}`, true); }
  finally { admin.submitting = false; renderFiles(); }
}

$('adminMode').addEventListener('change', () => {
  const pad = $('adminMode').value === 'four_side_pads';
  $('adminPadFields').hidden = !pad; $('adminOutletFields').hidden = pad;
});
$('refreshAdmin').addEventListener('click', refresh);
$('selectAll').addEventListener('click', () => {
  admin.selected = new Set(admin.files.filter(file => !active(file)).map(file => file.id));
  renderFiles();
});
$('clearAll').addEventListener('click', () => {
  admin.selected.clear(); renderFiles(); message('请选择至少一个可运行的 GDS。');
});
$('submitBatch').addEventListener('click', submit);
refresh();
setInterval(poll, 2000);
document.addEventListener('visibilitychange', () => { if (!document.hidden) poll(); });
