// The animations of the theory page (docs/theory.md): the track of a beta, the wall
// losses and the chamber seen from above. They are illustrations, apart from the
// wall-loss drawing, which repeats the calculation of ion_chambers.beta for a flat wall.
(function () {
if (!document.getElementById("ic-track")) return;

// colours and font shared by the three drawings, taken from the theme of the page
var pageColours = {};
var pageFont = getComputedStyle(document.body).fontFamily;
function readPageColours() {
  var style = getComputedStyle(document.body);
  ["mol", "beta", "pos", "neg", "exc", "panel", "muted", "ink", "line"].forEach(function (name) {
    pageColours[name] = style.getPropertyValue("--ic-" + name).trim();
  });
}
readPageColours();
var redrawHooks = [];
function onThemeChange() { readPageColours(); redrawHooks.forEach(function (hook) { hook(); }); }
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", onThemeChange);
new MutationObserver(onThemeChange).observe(document.body, {
  attributes: true, attributeFilter: ["data-theme"]
});
function fitCanvas(canvas) {
  var ratio = window.devicePixelRatio || 1;
  var context = canvas.getContext("2d");
  canvas.width = Math.round(canvas.clientWidth * ratio);
  canvas.height = Math.round(canvas.clientHeight * ratio);
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  return context;
}
var prefersReducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

// the same wall-loss calculation as the chamber model, for a flat wall
function createWallModel() {
  var Q = 18591, REST = 510999, ALPHA = 1 / 137.036;
  var i;

  function massRange(energy) {           // Katz-Penfold practical range, kg/m2
    var e = energy / 1e6;
    return 412 * Math.pow(e, 1.265 - 0.0954 * Math.log(e)) * 1e-2;
  }
  // tables to turn a remaining range back into an energy
  var tableEnergy = [], tableRange = [];
  for (i = 0; i < 400; i++) {
    var energy = Math.pow(10, i / 399 * Math.log10(Q));
    tableEnergy.push(energy); tableRange.push(massRange(energy));
  }
  function energyFromRange(range) {
    if (range <= tableRange[0]) return 0;
    var low = 0, high = tableRange.length - 1;
    if (range >= tableRange[high]) return tableEnergy[high];
    while (high - low > 1) {
      var middle = (low + high) >> 1;
      if (tableRange[middle] <= range) low = middle; else high = middle;
    }
    var share = (range - tableRange[low]) / (tableRange[high] - tableRange[low]);
    return tableEnergy[low] + share * (tableEnergy[high] - tableEnergy[low]);
  }
  // the tritium beta spectrum, as a cumulative distribution
  var gridEnergy = [], cumulative = [], total = 0;
  for (i = 0; i < 600; i++) {
    var E = 1 + i / 599 * (Q - 2);
    var W = E + REST, p = Math.sqrt(W * W - REST * REST);
    var eta = 2 * ALPHA * W / p;
    var fermi = 2 * Math.PI * eta / (1 - Math.exp(-2 * Math.PI * eta));
    total += fermi * p * W * (Q - E) * (Q - E);
    gridEnergy.push(E); cumulative.push(total);
  }
  // a fixed set of betas: energies from the spectrum, directions uniform on the sphere
  var seed = 12345;
  function random() {
    seed |= 0; seed = seed + 0x6D2B79F5 | 0;
    var t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  }
  var betas = [];
  for (i = 0; i < 3000; i++) {
    var target = random() * total, index = 0;
    while (cumulative[index] < target) index++;
    var cosPolar = 2 * random() - 1, azimuth = 2 * Math.PI * random();
    var sinPolar = Math.sqrt(1 - cosPolar * cosPolar);
    betas.push({ energy: gridEnergy[index], massRange: massRange(gridEnergy[index]),
                 ux: sinPolar * Math.cos(azimuth), uy: sinPolar * Math.sin(azimuth) });
  }

  // distance to the wall in metres, gas density in kg/m3
  function evaluate(distance, density) {
    var tracks = [], deposited = 0, emitted = 0, hits = 0;
    for (var k = 0; k < betas.length; k++) {
      var beta = betas[k];
      var range = beta.massRange / density;
      var toWall = beta.ux > 0 ? distance / beta.ux : Infinity;
      var reaches = toWall < range;
      var length = reaches ? toWall : range;
      var energyAtWall = reaches ? energyFromRange((range - toWall) * density) : 0;
      deposited += beta.energy - energyAtWall;
      emitted += beta.energy;
      if (reaches) hits++;
      if (k < 220) tracks.push({ dx: length * beta.ux, dy: length * beta.uy, reaches: reaches });
    }
    return { tracks: tracks, efficiency: deposited / emitted, hitShare: hits / betas.length,
             longest: massRange(Q) / density, typical: massRange(5700) / density };
  }
  return { evaluate: evaluate };
}

// core of the chamber animation, shared with the page (copied there verbatim)
function createChamber(random) {
  var ROD = 0.03;          // rod radius / can radius, to scale
  var REACH = 0.04;        // how close two ions must be to recombine
  var chamber = { positive: [], negative: [], flashes: [], created: 0, collected: 0,
                  recombined: 0, pending: 0, nextId: 1, ROD: ROD };

  chamber.reset = function () {
    chamber.positive = []; chamber.negative = []; chamber.flashes = [];
    chamber.created = 0; chamber.collected = 0; chamber.recombined = 0; chamber.pending = 0;
  };

  chamber.step = function (dt, pairsPerSecond, voltage, recombinationRate) {
    // 1. new ion pairs, uniformly over the gas
    chamber.pending += pairsPerSecond * dt;
    while (chamber.pending >= 1) {
      chamber.pending -= 1;
      var r = Math.sqrt(ROD * ROD + random() * (1 - ROD * ROD));
      var angle = random() * 2 * Math.PI;
      var id = chamber.nextId++;
      chamber.positive.push({ r: r, angle: angle, id: id });
      chamber.negative.push({ r: r, angle: angle, id: id });
      chamber.created += 1;
    }
    // 2. drift: the speed follows the field, which goes as 1/r
    var drift = 0.1 * voltage / 500;
    var i, ion;
    for (i = 0; i < chamber.positive.length; i++) {
      ion = chamber.positive[i];
      ion.r += drift / ion.r * dt;
      ion.angle += (random() - 0.5) * 0.25 * Math.sqrt(dt);
    }
    for (i = 0; i < chamber.negative.length; i++) {
      ion = chamber.negative[i];
      ion.r -= 1.1 * drift / ion.r * dt;
      ion.angle += (random() - 0.5) * 0.25 * Math.sqrt(dt);
    }
    // 3. collection at the electrodes
    var before = chamber.positive.length;
    chamber.positive = chamber.positive.filter(function (p) { return p.r < 1; });
    chamber.collected += before - chamber.positive.length;
    chamber.negative = chamber.negative.filter(function (n) { return n.r > ROD; });
    // 4. recombination of a positive and a negative ion that are close
    var cells = {};
    for (i = 0; i < chamber.negative.length; i++) {
      ion = chamber.negative[i];
      ion.x = ion.r * Math.cos(ion.angle); ion.y = ion.r * Math.sin(ion.angle);
      var key = Math.floor(ion.x / REACH) + "," + Math.floor(ion.y / REACH);
      (cells[key] = cells[key] || []).push(ion);
    }
    var chance = recombinationRate * dt;
    for (i = 0; i < chamber.positive.length; i++) {
      var p = chamber.positive[i];
      var px = p.r * Math.cos(p.angle), py = p.r * Math.sin(p.angle);
      var cx = Math.floor(px / REACH), cy = Math.floor(py / REACH);
      for (var a = -1; a <= 1 && !p.gone; a++) {
        for (var b = -1; b <= 1 && !p.gone; b++) {
          var list = cells[(cx + a) + "," + (cy + b)];
          if (!list) continue;
          for (var k = 0; k < list.length; k++) {
            var n = list[k];
            if (n.gone || n.id === p.id) continue;
            if (Math.hypot(n.x - px, n.y - py) < REACH && random() < chance) {
              n.gone = true; p.gone = true;
              chamber.recombined += 1;
              chamber.flashes.push({ x: px, y: py, age: 0 });
              break;
            }
          }
        }
      }
    }
    chamber.positive = chamber.positive.filter(function (q) { return !q.gone; });
    chamber.negative = chamber.negative.filter(function (q) { return !q.gone; });
    for (i = 0; i < chamber.flashes.length; i++) chamber.flashes[i].age += dt;
    chamber.flashes = chamber.flashes.filter(function (f) { return f.age < 0.5; });
  };
  return chamber;
}

(function () {
  var canvas = document.getElementById("ic-track");
  var ctx = canvas.getContext("2d");
  var energyEl = document.getElementById("ic-energy");
  var meterEl = document.getElementById("ic-meter");
  var pairsEl = document.getElementById("ic-pairs");
  var perPairEl = document.getElementById("ic-per-pair");
  var statusEl = document.getElementById("ic-status");
  var playButton = document.getElementById("ic-play");

  var START_ENERGY = 5700;      // eV, the mean energy of a tritium beta
  var IONISATION_ENERGY = 15.6; // eV, to free an electron from N2
  var EXCITATION_ENERGY = 10;   // eV, given to a molecule that is only excited
  var MEAN_ELECTRON_ENERGY = 21; // eV, mean kinetic energy of a freed electron
  var SECOND_PAIR_THRESHOLD = 34; // eV, above this a freed electron makes one more pair
  var SPACING = 46;             // px between molecules
  var REACH = 30;               // px, how close the beta must pass to interact

  var width = 0, height = 0;
  var state, speed = 4, running = true, lastTime = null;
  
  function resize() {
    var ratio = window.devicePixelRatio || 1;
    width = canvas.clientWidth;
    height = canvas.clientHeight;
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  }

  function newState() {
    state = {
      energy: START_ENERGY, pairs: 0, time: 0, worldX: 0, generatedTo: 0,
      betaY: height / 2, betaVy: 0, molecules: [], electrons: [], trail: [], done: false
    };
    addMolecules();
  }

  function addMolecules() {
    while (state.generatedTo < state.worldX + width + SPACING) {
      for (var y = SPACING / 2; y < height; y += SPACING) {
        if (Math.random() < 0.12) continue;
        state.molecules.push({
          x: state.generatedTo + (Math.random() - 0.5) * 24,
          y: y + (Math.random() - 0.5) * 24,
          angle: Math.random() * Math.PI,
          kind: "neutral", passed: false, excitedAt: -10
        });
      }
      state.generatedTo += SPACING;
    }
  }

  function freeElectron(molecule, kineticEnergy, awayFromY) {
    var direction = molecule.y >= awayFromY ? 1 : -1;
    var angle = direction * (0.5 + Math.random() * 0.9);
    var launch = 60 + 14 * Math.sqrt(kineticEnergy);
    var electron = {
      x: molecule.x, y: molecule.y,
      vx: launch * Math.cos(angle) * 0.6, vy: launch * Math.sin(angle),
      target: null
    };
    if (kineticEnergy > SECOND_PAIR_THRESHOLD) {
      // fast enough to ionise one more molecule: head for the nearest neutral one
      var best = null, bestDistance = 120;
      state.molecules.forEach(function (other) {
        if (other.kind !== "neutral" || other === molecule) return;
        var distance = Math.hypot(other.x - molecule.x, other.y - molecule.y);
        if (distance < bestDistance) { best = other; bestDistance = distance; }
      });
      if (best) {
        electron.target = best;
        electron.vx = (best.x - molecule.x) / bestDistance * launch;
        electron.vy = (best.y - molecule.y) / bestDistance * launch;
      }
    }
    state.electrons.push(electron);
  }

  function meet(molecule) {
    if (Math.abs(molecule.y - state.betaY) > REACH) return;
    var roll = Math.random();
    if (roll < 0.5 && molecule.kind === "neutral") {
      var room = state.energy - IONISATION_ENERGY;
      var kineticEnergy = Math.min(-MEAN_ELECTRON_ENERGY * Math.log(1 - Math.random()), room);
      state.energy -= IONISATION_ENERGY + kineticEnergy;
      state.pairs += 1;
      molecule.kind = "ion";
      freeElectron(molecule, kineticEnergy, state.betaY);
      state.betaVy += (Math.random() - 0.5) * 50;
    } else if (roll < 0.7 && molecule.kind === "neutral") {
      state.energy -= Math.min(EXCITATION_ENERGY, state.energy);
      molecule.excitedAt = state.time;
    }
    if (state.energy < IONISATION_ENERGY) {
      state.energy = 0;
      state.done = true;
    }
  }

  function step(dt) {
    if (state.done) { moveElectrons(dt); state.time += dt; return; }
    state.time += dt;
    var betaSpeed = 150 * (0.3 + 0.7 * Math.sqrt(state.energy / START_ENERGY));
    state.worldX += betaSpeed * dt;
    state.betaVy += (height / 2 - state.betaY) * 1.5 * dt;
    state.betaVy *= Math.pow(0.3, dt);
    state.betaY += state.betaVy * dt;
    state.betaY = Math.max(height * 0.25, Math.min(height * 0.75, state.betaY));

    var betaX = state.worldX + width * 0.3;
    state.trail.push({ x: betaX, y: state.betaY });
    if (state.trail.length > 60) state.trail.shift();

    addMolecules();
    for (var i = 0; i < state.molecules.length && !state.done; i++) {
      var molecule = state.molecules[i];
      if (!molecule.passed && molecule.x <= betaX) {
        molecule.passed = true;
        meet(molecule);
      }
    }
    state.molecules = state.molecules.filter(function (m) { return m.x > state.worldX - 60; });
    moveElectrons(dt);
  }

  function moveElectrons(dt) {
    state.electrons.forEach(function (electron) {
      electron.x += electron.vx * dt;
      electron.y += electron.vy * dt;
      var target = electron.target;
      if (target && Math.hypot(target.x - electron.x, target.y - electron.y) < 10) {
        electron.target = null;
        if (target.kind === "neutral") {
          target.kind = "ion";
          state.pairs += 1;
          freeElectron(target, 5, electron.y);
        }
      }
      if (!electron.target) {
        var slowing = Math.pow(0.04, dt);
        electron.vx *= slowing;
        electron.vy *= slowing;
      }
    });
    state.electrons = state.electrons.filter(function (e) { return e.x > state.worldX - 60; });
  }

  function drawMolecule(molecule, x) {
    var dx = Math.cos(molecule.angle) * 5, dy = Math.sin(molecule.angle) * 5;
    ctx.fillStyle = molecule.kind === "ion" ? pageColours.pos : pageColours.mol;
    ctx.beginPath();
    ctx.arc(x - dx, molecule.y - dy, 7, 0, 2 * Math.PI);
    ctx.arc(x + dx, molecule.y + dy, 7, 0, 2 * Math.PI);
    ctx.fill();
    if (molecule.kind === "ion") {
      ctx.fillStyle = pageColours.panel;
      ctx.font = "600 13px " + pageFont;
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("+", x, molecule.y + 0.5);
    }
    var sinceExcited = state.time - molecule.excitedAt;
    if (sinceExcited < 0.9) {
      ctx.globalAlpha = 1 - sinceExcited / 0.9;
      ctx.strokeStyle = pageColours.exc;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.arc(x, molecule.y, 14 + sinceExcited * 10, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.globalAlpha = 1;
    }
  }

  function draw() {
    ctx.clearRect(0, 0, width, height);
    var offset = state.worldX;

    state.molecules.forEach(function (molecule) { drawMolecule(molecule, molecule.x - offset); });

    ctx.fillStyle = pageColours.neg;
    state.electrons.forEach(function (electron) {
      ctx.beginPath();
      ctx.arc(electron.x - offset, electron.y, 3.2, 0, 2 * Math.PI);
      ctx.fill();
    });

    ctx.strokeStyle = pageColours.beta;
    ctx.lineWidth = 2;
    ctx.setLineDash([4, 5]);
    ctx.beginPath();
    state.trail.forEach(function (point, i) {
      if (i === 0) ctx.moveTo(point.x - offset, point.y);
      else ctx.lineTo(point.x - offset, point.y);
    });
    ctx.stroke();
    ctx.setLineDash([]);

    var last = state.trail[state.trail.length - 1];
    if (last) {
      ctx.fillStyle = pageColours.beta;
      ctx.shadowColor = pageColours.beta;
      ctx.shadowBlur = state.done ? 0 : 14;
      ctx.beginPath();
      ctx.arc(last.x - offset, last.y, 6, 0, 2 * Math.PI);
      ctx.fill();
      ctx.shadowBlur = 0;
    }
  }

  function thousands(value) {
    return String(Math.round(value)).replace(/\B(?=(\d{3})+(?!\d))/g, " ");
  }

  function updateReadouts() {
    energyEl.textContent = thousands(state.energy) + " eV";
    meterEl.style.width = (100 * state.energy / START_ENERGY) + "%";
    pairsEl.textContent = state.pairs;
    var spent = START_ENERGY - state.energy;
    perPairEl.textContent = state.pairs > 4 ? (spent / state.pairs).toFixed(1) + " eV" : "–";
    if (state.done) {
      statusEl.textContent = "The beta has stopped: " + state.pairs + " ion pairs from 5 700 eV.";
      playButton.disabled = true;
    } else {
      statusEl.textContent = running ? "" : "Paused.";
      playButton.disabled = false;
    }
    playButton.textContent = running && !state.done ? "Pause" : "Play";
  }

  function frame(now) {
    if (lastTime === null) lastTime = now;
    var elapsed = Math.min((now - lastTime) / 1000, 0.05);
    lastTime = now;
    if (running) {
      var steps = speed;
      for (var i = 0; i < steps; i++) step(elapsed);
    }
    draw();
    updateReadouts();
    requestAnimationFrame(frame);
  }

  function start() {
    newState();
    // show a few events on the first frame
    for (var i = 0; i < 150; i++) step(1 / 60);
  }

  playButton.addEventListener("click", function () {
    running = !running;
  });
  document.getElementById("ic-restart").addEventListener("click", function () {
    start();
    running = true;
  });
  [1, 4, 16].forEach(function (value) {
    document.getElementById("ic-speed" + value).addEventListener("click", function () {
      speed = value;
      [1, 4, 16].forEach(function (other) {
        document.getElementById("ic-speed" + other).setAttribute("aria-pressed", String(other === value));
      });
    });
  });

  window.addEventListener("resize", function () { resize(); });

  resize();
  start();
  if (prefersReducedMotion) running = false;
  requestAnimationFrame(frame);
})();

(function () {
  var model = createWallModel();
  var canvas = document.getElementById("ic-wall-canvas");
  var distanceInput = document.getElementById("ic-wall-distance");
  var pressureInput = document.getElementById("ic-wall-pressure");
  var AIR_DENSITY = 1.205; // kg/m3 at 1 atm

  function draw() {
    var context = fitCanvas(canvas);
    var width = canvas.clientWidth, height = canvas.clientHeight;
    var distance = Number(distanceInput.value) * 1e-3;
    var pressure = Number(pressureInput.value);
    var result = model.evaluate(distance, AIR_DENSITY * pressure);

    document.getElementById("ic-wall-distance-out").textContent = (distance * 1e3).toFixed(2) + " mm";
    document.getElementById("ic-wall-pressure-out").textContent = pressure.toFixed(1) + " atm";
    document.getElementById("ic-wall-hits").textContent = Math.round(100 * result.hitShare) + " %";
    document.getElementById("ic-wall-kept").textContent = (100 * result.efficiency).toFixed(1) + " %";
    document.getElementById("ic-wall-range").textContent = (result.longest * 1e3).toFixed(1) + " mm";

    var wallX = width - 46;
    var scale = (wallX - 16) / 6e-3;            // pixels per metre: 6 mm of gas shown
    var originX = wallX - distance * scale, originY = height / 2;

    context.clearRect(0, 0, width, height);
    // the wall
    context.fillStyle = pageColours.mol;
    context.fillRect(wallX, 0, width - wallX, height);
    context.fillStyle = pageColours.ink;
    context.font = "500 12px " + pageFont;
    context.textAlign = "center";
    context.save();
    context.translate(wallX + 24, height / 2);
    context.rotate(-Math.PI / 2);
    context.fillText("metal wall", 0, 0);
    context.restore();
    // a 1 mm scale bar
    context.strokeStyle = pageColours.muted;
    context.lineWidth = 2;
    context.beginPath();
    context.moveTo(16, height - 18);
    context.lineTo(16 + 1e-3 * scale, height - 18);
    context.stroke();
    context.fillStyle = pageColours.muted;
    context.textAlign = "left";
    context.fillText("1 mm", 16, height - 26);

    // the tracks, clipped to the gas
    context.save();
    context.beginPath();
    context.rect(0, 0, wallX, height);
    context.clip();
    context.lineWidth = 1;
    for (var pass = 0; pass < 2; pass++) {
      context.strokeStyle = pass === 0 ? pageColours.muted : pageColours.pos;
      context.globalAlpha = pass === 0 ? 0.55 : 0.9;
      context.beginPath();
      result.tracks.forEach(function (track) {
        if (track.reaches !== (pass === 1)) return;
        context.moveTo(originX, originY);
        context.lineTo(originX + track.dx * scale, originY + track.dy * scale);
      });
      context.stroke();
    }
    context.restore();
    context.globalAlpha = 1;
    context.fillStyle = pageColours.beta;
    context.beginPath();
    context.arc(originX, originY, 5, 0, 2 * Math.PI);
    context.fill();
  }

  distanceInput.addEventListener("input", draw);
  pressureInput.addEventListener("input", draw);
  window.addEventListener("resize", draw);
  redrawHooks.push(draw);
  draw();
})();

(function () {
  var canvas = document.getElementById("ic-chamber-canvas");
  var chamber = createChamber(Math.random);
  var playButton = document.getElementById("ic-sim-play");
  var RECOMBINATION_RATE = 0.3;
  var pairsPerSecond = 60, voltage = 500, running = !prefersReducedMotion, lastTime = null;
  var context = fitCanvas(canvas);

  function restart() {
    chamber.reset();
    // start with ions already crossing the gas, then count from zero
    for (var i = 0; i < 360; i++) chamber.step(1 / 60, pairsPerSecond, voltage, RECOMBINATION_RATE);
    chamber.created = 0; chamber.collected = 0; chamber.recombined = 0;
  }

  function draw() {
    var width = canvas.clientWidth, height = canvas.clientHeight;
    var radius = Math.min(width, height) / 2 - 22;
    var centreX = width / 2, centreY = height / 2;
    context.clearRect(0, 0, width, height);

    context.strokeStyle = pageColours.ink;
    context.lineWidth = 3;
    context.beginPath();
    context.arc(centreX, centreY, radius, 0, 2 * Math.PI);
    context.stroke();

    function drawIons(ions, colour) {
      context.fillStyle = colour;
      context.beginPath();
      ions.forEach(function (ion) {
        var x = centreX + radius * ion.r * Math.cos(ion.angle);
        var y = centreY + radius * ion.r * Math.sin(ion.angle);
        context.moveTo(x + 2.6, y);
        context.arc(x, y, 2.6, 0, 2 * Math.PI);
      });
      context.fill();
    }
    drawIons(chamber.positive, pageColours.pos);
    drawIons(chamber.negative, pageColours.neg);

    context.strokeStyle = pageColours.exc;
    context.lineWidth = 2;
    chamber.flashes.forEach(function (flash) {
      context.globalAlpha = 1 - flash.age / 0.5;
      context.beginPath();
      context.arc(centreX + radius * flash.x, centreY + radius * flash.y, 4 + 22 * flash.age, 0, 2 * Math.PI);
      context.stroke();
    });
    context.globalAlpha = 1;

    context.fillStyle = pageColours.ink;
    context.beginPath();
    context.arc(centreX, centreY, Math.max(radius * chamber.ROD, 3), 0, 2 * Math.PI);
    context.fill();

    context.font = "500 12px " + pageFont;
    context.textAlign = "left";
    context.fillStyle = pageColours.muted;
    context.fillText("can, 0 V", 12, 20);
    context.fillText("rod, +" + voltage + " V", centreX + 10, centreY - 8);
  }

  function updateReadouts() {
    document.getElementById("ic-sim-collected").textContent = chamber.collected;
    document.getElementById("ic-sim-recombined").textContent = chamber.recombined;
    var finished = chamber.collected + chamber.recombined;
    document.getElementById("ic-sim-efficiency").textContent =
      finished >= 20 ? Math.round(100 * chamber.collected / finished) + " %" : "–";
    playButton.textContent = running ? "Pause" : "Play";
  }

  function frame(now) {
    if (lastTime === null) lastTime = now;
    var elapsed = Math.min((now - lastTime) / 1000, 0.05);
    lastTime = now;
    if (running) chamber.step(elapsed, pairsPerSecond, voltage, RECOMBINATION_RATE);
    draw();
    updateReadouts();
    requestAnimationFrame(frame);
  }

  function choose(prefix, values, chosen) {
    values.forEach(function (value) {
      document.getElementById(prefix + value).setAttribute("aria-pressed", String(value === chosen));
    });
  }
  [15, 60, 240].forEach(function (value) {
    document.getElementById("ic-act" + value).addEventListener("click", function () {
      pairsPerSecond = value; choose("ic-act", [15, 60, 240], value); restart();
    });
  });
  [200, 500, 1500].forEach(function (value) {
    document.getElementById("ic-volt" + value).addEventListener("click", function () {
      voltage = value; choose("ic-volt", [200, 500, 1500], value); restart();
    });
  });
  playButton.addEventListener("click", function () { running = !running; });
  window.addEventListener("resize", function () { context = fitCanvas(canvas); });

  restart();
  requestAnimationFrame(frame);
})();
})();
