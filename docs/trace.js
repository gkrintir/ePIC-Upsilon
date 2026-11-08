(() => {
  'use strict';

  function redirectFigure() {
    if (!/^#(?:mass|efficiency|w|differential|figure-content)(?:\?.*)?$/.test(location.hash)) return false;
    location.replace('./figures.html' + location.search + location.hash);
    return true;
  }
  if (redirectFigure()) return;

  const data = window.UPSiLON_TRACE;
  const byId = name => document.getElementById(`trace-${name}`);
  const status = byId('status');
  if (!data || !data.routes || !data.sources || !Object.keys(data.routes).length) {
    status.prepend('The interactive guide could not load. ');
    return;
  }

  const routeKeys = Object.keys(data.routes);
  const defaultRoute = routeKeys.includes('mass') ? 'mass' : routeKeys[0];
  const select = byId('result');
  for (const key of routeKeys) select.add(new Option(data.routes[key].label, key));

  function locationState() {
    const [key, stepId, mode] = location.hash.slice(1).split('?')[0].split('/');
    const result = routeKeys.includes(key) ? key : defaultRoute;
    const steps = data.routes[result].steps;
    // Older review links without a stage start at the published result.
    const step = steps.find(item => item.id === stepId) || (mode === 'review' ? steps.at(-1) : steps[0]);
    return {result, step};
  }

  function hash(result, stepId) {
    const params=figureParams(result);
    return `#${result}/${stepId}${params ? '?'+params : ''}`;
  }

  function figureParams(key) {
    const [route, query] = location.hash.slice(1).split("?");
    if (route.split("/")[0]!==key || !query) return "";
    const allowed = new Set(["state","group","observable","map","bin","x","y","scale","panel"]);
    const params = new URLSearchParams();
    for (const [name,value] of new URLSearchParams(query)) {
      if (allowed.has(name) && /^[a-zA-Z0-9_-]{1,55}$/.test(value)) params.set(name,value);
    }
    return params.toString();
  }
  function resultURL(route,key) {
    const path=data.sources[route.figure].url, params=figureParams(key);
    if (!params) return path;
    return path+(key==="beams" ? "#" : "?")+params;
  }

  function link(label, url) {
    const anchor = document.createElement('a');
    anchor.textContent = label;
    anchor.href = url;
    return anchor;
  }

  function render() {
    if (redirectFigure()) return;
    const {result, step} = locationState();
    const route = data.routes[result];
    const steps = route.steps;
    const position = steps.indexOf(step);
    const atStage = item => hash(result, item.id);
    select.value = result;
    byId('route-title').textContent = route.label;
    byId('description').textContent = route.description;
    const figureURL=resultURL(route,result);
    byId("references").replaceChildren(link(route.section+" · "+route.paper,figureURL));
    const boundary=byId("boundary");
    boundary.hidden=!["reconstruction","tuple"].includes(step.id);
    boundary.textContent=data.boundary;
    const shortcut=byId("shortcut");
    shortcut.replaceChildren();
    if (step.id==="reconstruction") {
      shortcut.append(link("Already have merged ROOT files? Start here →",atStage(steps[1])));
    } else if (step.id==="tuple" || step.id==="analysis") {
      const run=data.sources[result==="beams" ? "beamGuide" : "paperRun"];
      shortcut.append(link("Run this analysis →",run.url));
    }
    const preview=byId("preview");
    preview.replaceChildren();
    if (step.id==="publication" && route.preview) {
      const thumb=link("",figureURL), image=document.createElement("img");
      image.src=route.preview; image.alt="Open "+route.paper+" in the figure viewer";thumb.append(image);
      preview.append(thumb);
    }


    const stageLinks = steps.map((item, index) => {
      const li = document.createElement('li');
      const anchor = link(item.title, atStage(item));
      if (item.id === step.id) anchor.setAttribute('aria-current', 'step');
      const badge = document.createElement('span');
      badge.className = 'trace-step-number';
      badge.setAttribute('aria-hidden', 'true');
      badge.textContent = String(index + 1).padStart(2, '0');
      anchor.prepend(badge);
      li.append(anchor);
      return li;
    });
    byId('stages').replaceChildren(...stageLinks);
    const progress = `Stage ${position + 1} of ${steps.length}`;
    byId('progress').textContent = progress;
    byId('title').textContent = step.title;
    byId('detail').textContent = step.detail;
    byId('input').textContent = step.input;
    byId('output').textContent = step.output;
    const sources = step.links.map(key => data.sources[key]).filter(Boolean).map(source => {
      const li = document.createElement('li');
      li.append(link(source.label,source === data.sources[route.figure] ? figureURL : source.url));
      return li;
    });
    byId('sources').replaceChildren(...sources);
    const neighbours = [];
    if (position > 0) {
      const prev = steps[position - 1];
      neighbours.push(link(`← ${prev.title}`, atStage(prev)));
    }
    if (position < steps.length - 1) {
      const next = steps[position + 1];
      const nextLink = link(`${next.title} →`, atStage(next));
      nextLink.className = 'trace-next';
      neighbours.push(nextLink);
    }
    byId('pagination').replaceChildren(...neighbours);
    byId('announcer').textContent = `${route.label}. ${progress}. ${step.title}.`;
    window.showTraceRun?.(result, step.id);
    document.title = `${step.title} · ${route.label} · ePIC Υ Atlas`;
  }

  select.addEventListener('change', () => {
    location.hash = hash(select.value, locationState().step.id);
  });
  let focusStage = false;
  for (const id of ['stages', 'pagination']) {
    byId(id).addEventListener('click', event => {
      const anchor = event.target.closest('a');
      if (!anchor || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
      if (anchor.hash === location.hash) byId('title').focus();
      else focusStage = true;
    });
  }
  // Keep the skip link independent of the URL fragment that identifies the current stage.
  document.querySelector('.skip').addEventListener('click', event => {
    event.preventDefault();
    byId('main').focus();
  });
  window.addEventListener('hashchange', () => {
    render();
    if (focusStage) byId('title').focus();
    focusStage = false;
  });
  render();
  status.hidden = true;
  byId('content').hidden = false;
})();
