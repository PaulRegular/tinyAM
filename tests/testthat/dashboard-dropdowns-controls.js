// A small DOM/Bootstrap contract double: validates template adaptation, not browser layout.
const fs = require('fs'), assert = require('assert');
class Element {
  constructor(tag) {
    this.tagName = tag; this.children = []; this.dataset = {}; this.attributes = {};
    this.events = {}; this.classes = new Set(); this.textContent = '';
    this.offsetWidth = 100; this.offsetHeight = 100;
    this.classList = {contains: c => this.classes.has(c), toggle: (c, yes) =>
      yes ? this.classes.add(c) : this.classes.delete(c)};
  }
  set className(value) { this.classes = new Set(value.split(' ')); }
  get className() { return Array.from(this.classes).join(' '); }
  appendChild(child) { child.parentElement = this; this.children.push(child); }
  insertBefore(child, target) {
    child.parentElement = this; this.children.splice(this.children.indexOf(target), 0, child);
  }
  setAttribute(key, value) { this.attributes[key] = value; }
  getAttribute(key) { return this.attributes[key]; }
  addEventListener(key, fn) { this.events[key] = fn; }
  focus() { document.activeElement = this; }
  all() { return this.children.flatMap(child => [child, ...child.all()]); }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  querySelectorAll(selector) {
    if (selector === ':scope > .nav-tabs') return this.children.filter(x => x.classes.has('nav-tabs'));
    if (selector === ':scope > .tab-content > [data-tab-group]')
      return this.children.filter(x => x.classes.has('tab-content')).flatMap(x => x.children).filter(x => x.dataset.tabGroup);
    return this.all().filter(x => selector === '[data-dropdown-tabs]' ? x.dataset.dropdownTabs :
      selector === 'a[role="tab"]' ? x.tagName === 'a' && x.attributes.role === 'tab' :
      selector === 'button' ? x.tagName === 'button' :
      selector === '.js-plotly-plot' ? x.classes.has('js-plotly-plot') : false);
  }
}
global.document = new Element('document');
document.createElement = tag => new Element(tag);
document.getElementById = id => document.all().find(x => x.id === id);
global.window = global;
let resizes = 0;
global.Plotly = {Plots:{resize:() => { resizes++; }}};
global.requestAnimationFrame = fn => fn();
global.jQuery = target => {
  if (typeof target === 'function') { target(); return; }
  return {on:(event, fn) => { target.events[event] = fn; },
    dropdown:() => { target.setAttribute('aria-expanded', 'false'); },
    tab:() => {
      const nav = target.parentElement.parentElement;
      nav.querySelectorAll('a[role="tab"]').forEach(a => a.classList.toggle('active', a === target));
      nav.events['shown.bs.tab']({target});
    }};
};
function row(prefix, initial) {
  const row = new Element('div'); row.dataset.dropdownTabs = 'true';
  const nav = new Element('ul'); nav.className = 'nav nav-tabs';
  const content = new Element('div'); content.className = 'tab-content';
  row.appendChild(nav); row.appendChild(content); document.appendChild(row);
  const charts = [['Recruitment'], ...(initial ? [['States','N'],['Initial abundance','N']] : [['N']]),
    ['F'], ['M'], ['Recruitment','Process deviations'], ['N/F/M','Process deviations'],
    ['Recruitment','Joint-sample residuals'], ['N/F/M','Joint-sample residuals'], ['Formula effect']];
  charts.forEach(([label,group], i) => {
    const pane = new Element('div'); pane.id = prefix + i;
    if (group) pane.dataset.tabGroup = group;
    const plot = new Element('div'); plot.className = 'js-plotly-plot'; pane.appendChild(plot);
    content.appendChild(pane);
    const li = new Element('li'), anchor = new Element('a'); anchor.textContent = label;
    anchor.setAttribute('role','tab'); anchor.setAttribute('href','#'+pane.id);
    anchor.setAttribute('aria-controls',pane.id); anchor.classList.toggle('active',i===0);
    li.appendChild(anchor); nav.appendChild(li);
  });
  return {row,nav,content};
}
const a = row('a', true), b = row('b', false);
eval(fs.readFileSync(process.argv[2], 'utf8'));
function visible(nav) { return nav.children.filter(li => !li.hidden); }
assert.deepEqual(visible(a.nav).map(li => li.children[0].textContent),
  ['Recruitment','N','F','M','Process deviations','Joint-sample residuals','Formula effect']);
assert.equal(a.content.children.length, 10);
const n = visible(a.nav)[1], processMenu = visible(a.nav)[4], residual = visible(a.nav)[5];
assert.equal(processMenu.children[0].tagName, 'a');
assert.equal(processMenu.children[0].attributes.role, 'button');
assert.equal(processMenu.children[0].getAttribute('href'), '#');
assert.deepEqual(n.children[1].children.map(x => x.textContent), ['States','Initial abundance']);
assert.deepEqual(processMenu.children[1].children.map(x => x.textContent), ['Recruitment','N/F/M']);
assert.deepEqual(residual.children[1].children.map(x => x.textContent), ['Recruitment','N/F/M']);
assert.equal(visible(b.nav)[1].children.length, 1); // N alone has no dropdown
processMenu.children[1].children[1].events.click();
assert(processMenu.children[0].classList.contains('dropdown-selected'));
assert(processMenu.children[1].children[1].classList.contains('dropdown-selected'));
assert.equal(resizes, 1);
assert(!visible(b.nav)[4].children[0].classList.contains('dropdown-selected'));
residual.children[1].children[0].events.click();
assert(!processMenu.children[0].classList.contains('dropdown-selected'));
assert(processMenu.children[1].children[1].classList.contains('dropdown-selected'));
assert(residual.children[0].classList.contains('dropdown-selected'));
const menu = n.children[1]; menu.children[0].focus();
menu.events.keydown({key:'ArrowDown',preventDefault(){},stopPropagation(){}});
assert.equal(document.activeElement, menu.children[1]);
menu.events.keydown({key:'ArrowUp',preventDefault(){},stopPropagation(){}});
assert.equal(document.activeElement, menu.children[0]);
menu.events.keydown({key:'Escape',preventDefault(){},stopPropagation(){}});
assert.equal(document.activeElement,n.children[0]);
document.events['flexdashboard:layoutcomplete']();
assert.equal(visible(a.nav).length,7); // repeated layout does not duplicate menus
assert.equal(a.content.children.length,10);
console.log('Dashboard dropdown controls passed');
