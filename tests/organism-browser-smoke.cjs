// Optional live smoke test: serve the app on 8765 and start headless Chromium/Edge
// with --remote-debugging-port=9222 and an isolated temporary user-data-dir.
// Node 22+; no packages. Outputs screenshots/downloads to the OS temporary folder.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const base = process.env.ATLAS_URL || 'http://127.0.0.1:8765/';
const debug = process.env.ATLAS_CDP_URL || 'http://127.0.0.1:9222';
const output = fs.mkdtempSync(path.join(os.tmpdir(), 'aa-atlas-organisms-'));

(async () => {
  const tab = await (await fetch(`${debug}/json/new?about:blank`, {method:'PUT'})).json();
  const ws = new WebSocket(tab.webSocketDebuggerUrl);
  await new Promise((resolve,reject)=>{ws.onopen=resolve;ws.onerror=reject;});
  const pending = new Map(), exceptions = [], requests = [];
  let next = 0, missing = false;
  const send = (method, params={}) => new Promise((resolve,reject)=>{
    const id=++next;pending.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params}));
  });
  ws.onmessage = event => {
    const m=JSON.parse(event.data);
    if(m.id){const p=pending.get(m.id);pending.delete(m.id);m.error?p.reject(Error(JSON.stringify(m.error))):p.resolve(m.result);}
    if(m.method==='Runtime.exceptionThrown') exceptions.push(m.params.exceptionDetails);
    if(m.method==='Network.requestWillBeSent') requests.push(m.params.request.url);
    if(m.method==='Fetch.requestPaused') {
      const {requestId,request}=m.params;
      const optional=/annotations\/ecoli\/|Ec_results\/D\/manifest|Ec_results\/electrostatics\/|Ec_results\/pocket_points|Ec_results\/L\/aa_positions\//.test(request.url);
      void send(missing && optional?'Fetch.fulfillRequest':'Fetch.continueRequest',missing && optional?{requestId,responseCode:404,body:''}:{requestId});
    }
  };
  const evaluate = async expression => {
    const r=await send('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});
    if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));
    return r.result.value;
  };
  const wait = async (expression, timeout=120000) => {
    const start=Date.now();
    while(Date.now()-start<timeout){try{if(await evaluate(expression))return;}catch{}
      await new Promise(resolve=>setTimeout(resolve,100));}
    throw Error(`Timed out: ${expression}`);
  };
  const ready = organism => wait(`typeof ORGANISM!=='undefined' && ORGANISM.id==='${organism}' && !organismLeaving && document.querySelector('#loading-screen').classList.contains('hidden')`);
  const select = organism => evaluate(`document.querySelector('#organism-select').value='${organism}';document.querySelector('#organism-select').dispatchEvent(new Event('change'));`);
  const click = selector => evaluate(`document.querySelector(${JSON.stringify(selector)}).click()`);
  try {
    await send('Runtime.enable');await send('Network.enable');await send('Page.enable');
    await send('Page.setDownloadBehavior',{behavior:'allow',downloadPath:output});
    await send('Emulation.setDeviceMetricsOverride',{width:1440,height:1000,deviceScaleFactor:1,mobile:false});
    await send('Page.navigate',{url:base});await ready('ecoli');
    assert.equal(await evaluate('state.metadata.size'),3085);
    assert.match(await evaluate(`document.querySelector('.hero .eyebrow').textContent`),/Escherichia coli/);
    assert.match(await evaluate(`document.querySelector('footer').textContent`),/E. coli/);
    assert.equal(await evaluate('state.rawByAA.size'),39);
    assert.ok(requests.filter(url=>/results|annotations/.test(url)).every(url=>!url.includes('At_results/')&&!url.includes('annotations/arabidopsis/')));
    await wait(`document.querySelector('#explorer-qphi-status').hidden`);
    await click('#download-csv');
    await wait(`!document.querySelector('#download-csv').disabled`);
    // Downloads are asynchronous outside the page context.
    for(let i=0;i<100&&!fs.existsSync(path.join(output,'ecoli_aa_atlas_ala_best_pockets_filtered.csv'));i++)await new Promise(r=>setTimeout(r,100));
    const csv=fs.readFileSync(path.join(output,'ecoli_aa_atlas_ala_best_pockets_filtered.csv'),'utf8');
    assert.match(csv,/gene_id/);assert.doesNotMatch(csv,/tair_id/);
    await click('[data-gene-info]');
    await wait(`document.querySelector('#gene-description-body').getAttribute('aria-busy')==='false'`);
    assert.match(await evaluate(`document.querySelector('#gene-description-title').textContent`),/E. coli/);
    assert.doesNotMatch(await evaluate(`document.querySelector('#gene-description-body').textContent`),/TAIR|Araport/);
    await click('#gene-description-close');
    for(const view of ['matrix','compare','overlap','statistics','methods','go','control-qc']) {
      await click(`[data-view="${view}"]`);
      if(view==='go')await wait('goResult!==null');
      if(view==='control-qc')await wait('controlQCResult!==null');
      assert.equal(await evaluate('state.currentView'),view);
    }
    await click('[data-view="protein"]');
    await wait(`document.querySelector('#structure-ligand-status').dataset.state==='ready'`,60000);
    await wait(`document.querySelector('#structure-cloud-status').dataset.state==='ready'`,60000);
    await evaluate(`document.querySelector('#pocket-point-color').value='potential';document.querySelector('#pocket-point-color').dispatchEvent(new Event('change'))`);
    await wait(`pocketPotentialStatus.kind==='ready'`,60000);
    console.log('E. coli viewer:',await evaluate(`({ligand:document.querySelector('#structure-ligand-status').textContent,cloud:document.querySelector('#structure-cloud-status').textContent})`));
    await click('[data-view="explorer"]');
    const screenshot=await send('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(output,'ecoli-desktop.png'),Buffer.from(screenshot.data,'base64'));
    await send('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:false});
    assert.ok(await evaluate(`(()=>{const r=document.querySelector('#organism-select').getBoundingClientRect();return r.left>=0&&r.right<=innerWidth})()`));
    const mobile=await send('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(output,'ecoli-mobile.png'),Buffer.from(mobile.data,'base64'));
    await send('Emulation.setDeviceMetricsOverride',{width:1440,height:1000,deviceScaleFactor:1,mobile:false});

    await select('arabidopsis');await ready('arabidopsis');
    assert.equal(await evaluate('state.rawByAA.get("ALA").length'),33622);
    assert.equal(await evaluate('goResult===null && controlQCResult===null && ligandPositionCache.size===0 && pocketPointTableCache===null'),true);
    assert.equal(await evaluate('state.annotations.has("P00561")'),false);
    await select('ecoli');await ready('ecoli');
    await evaluate(`switchOrganism('arabidopsis');switchOrganism('ecoli');switchOrganism('arabidopsis')`);
    await wait(`typeof ORGANISM!=='undefined'&&ORGANISM.id==='arabidopsis'&&!organismLeaving`);
    await select('ecoli');await ready('ecoli');
    assert.equal(await evaluate(`state.rawByAA.size===39 && [...state.rawByAA.values()].every(rows=>rows.every(row=>row._compactDirectory.startsWith('Ec_results/')))`),true);
    assert.equal(await evaluate(`state.search===''&&state.selectedPocket===null&&goResult===null&&controlQCResult===null`),true);

    missing=true;await send('Fetch.enable',{patterns:[{urlPattern:base+'*'}]});
    await send('Page.navigate',{url:base+'?organism=ecoli'});await ready('ecoli');
    assert.equal(await evaluate('state.rawByAA.size'),20);
    assert.match(await evaluate(`document.querySelector('#dataset-availability').textContent`),/Unavailable.*UniProt.*DALA/);
    await wait(`document.querySelector('#explorer-qphi-status').dataset.state==='unavailable'`);
    await click('[data-view="go"]');await wait(`!document.querySelector('#go-error').hidden`);
    assert.match(await evaluate(`document.querySelector('#go-error-message').textContent`),/Unavailable/);
    await click('[data-view="control-qc"]');await wait(`!document.querySelector('#control-qc-retry').hidden`);
    assert.match(await evaluate(`document.querySelector('#control-qc-status').textContent`),/unavailable/i);
    await click('[data-view="explorer"]');await click('[data-gene-info]');
    await wait(`document.querySelector('#gene-description-body').getAttribute('aria-busy')==='false'`);
    assert.match(await evaluate(`document.querySelector('#gene-description-body').textContent`),/Unavailable/);
    assert.equal(exceptions.length,0,JSON.stringify(exceptions));
    console.log(`Browser smoke passed: all tabs, both directions, rapid switches, downloads, optional 404s and mobile selector. Artifacts: ${output}`);
  } finally { ws.close(); }
})().catch(error=>{console.error(error);process.exitCode=1;});
