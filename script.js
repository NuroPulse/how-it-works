const $ = (id) => document.getElementById(id);
const fmt = (x, d = 1) => x.toFixed(d).replace(".", ",");
let seed = 11;
const rnd = () => {
	seed = (seed * 16807) % 2147483647;
	return seed / 2147483647;
};
const gauss = () => {
	let g = 0;
	for (let j = 0; j < 6; j++) g += rnd();
	return (g - 3) * 1.4142;
};
const T = 16,
	FS0 = 100,
	N0 = T * FS0;
const jit = [],
	nz = [];
for (let i = 0; i < N0; i++) {
	jit.push(i === 0 ? 0 : (rnd() - 0.5) * 0.9);
	nz.push(gauss());
}

function biquad(type, fc, fs, Q) {
	const w0 = (2 * Math.PI * fc) / fs,
		c = Math.cos(w0),
		s = Math.sin(w0),
		a = s / (2 * Q);
	let b0, b1, b2;
	if (type === "lp") {
		b0 = (1 - c) / 2;
		b1 = 1 - c;
		b2 = (1 - c) / 2;
	} else {
		b0 = (1 + c) / 2;
		b1 = -(1 + c);
		b2 = (1 + c) / 2;
	}
	const a0 = 1 + a;
	return {
		b0: b0 / a0,
		b1: b1 / a0,
		b2: b2 / a0,
		a1: (-2 * c) / a0,
		a2: (1 - a) / a0,
	};
}
function run(x, f) {
	const y = new Array(x.length);
	let x1 = 0,
		x2 = 0,
		y1 = 0,
		y2 = 0;
	for (let i = 0; i < x.length; i++) {
		const v = f.b0 * x[i] + f.b1 * x1 + f.b2 * x2 - f.a1 * y1 - f.a2 * y2;
		y[i] = v;
		x2 = x1;
		x1 = x[i];
		y2 = y1;
		y1 = v;
	}
	return y;
}
function filtfilt(x, f) {
	const a = run(x, f);
	a.reverse();
	const b = run(a, f);
	b.reverse();
	return b;
}
function resample(t, v, fs) {
	const t0 = t[0],
		n = Math.floor((t[t.length - 1] - t0) * fs) + 1,
		out = new Array(n);
	let j = 0;
	for (let i = 0; i < n; i++) {
		const tt = t0 + i / fs;
		while (j < t.length - 2 && t[j + 1] < tt) j++;
		const sp = t[j + 1] - t[j];
		const a = sp > 0 ? Math.min(1, Math.max(0, (tt - t[j]) / sp)) : 0;
		out[i] = v[j] + a * (v[j + 1] - v[j]);
	}
	return { x: out, t0: t0 };
}

const cosT = [],
	sinT = [];
for (let k = 0; k <= 60; k++) {
	const c = new Float64Array(200),
		s = new Float64Array(200);
	for (let n = 0; n < 200; n++) {
		const a = (2 * Math.PI * k * n) / 200;
		c[n] = Math.cos(a);
		s[n] = Math.sin(a);
	}
	cosT.push(c);
	sinT.push(s);
}
const win = new Float64Array(200);
let W = 0;
for (let n = 0; n < 200; n++) {
	win[n] = 0.5 * (1 - Math.cos((2 * Math.PI * n) / 200));
	W += win[n];
}

function welch(x) {
	const segAmp = [],
		P = new Float64Array(61);
	let ns = 0;
	for (let s = 0; s + 200 <= x.length; s += 100) {
		const w = new Float64Array(200);
		for (let n = 0; n < 200; n++) w[n] = x[s + n] * win[n];
		const a = new Float64Array(61);
		for (let k = 1; k <= 60; k++) {
			let re = 0,
				im = 0;
			const c = cosT[k],
				sn = sinT[k];
			for (let n = 0; n < 200; n++) {
				re += w[n] * c[n];
				im -= w[n] * sn[n];
			}
			const p = re * re + im * im;
			P[k] += p;
			a[k] = (2 * Math.sqrt(p)) / W;
		}
		segAmp.push(a);
		ns++;
	}
	const A = new Float64Array(61);
	for (let k = 1; k <= 60; k++)
		A[k] = (2 * Math.sqrt(P[k] / Math.max(ns, 1))) / W;
	return { A: A, segAmp: segAmp, ns: ns };
}

function proc() {
	const aT = +$("c1").value,
		aD = +$("c2").value,
		aN = +$("c3").value,
		fT = +$("c5").value,
		gap = $("c4").checked;
	$("c1v").textContent = fmt(aT, 2);
	$("c2v").textContent = fmt(aD, 1);
	$("c3v").textContent = fmt(aN, 2);
	$("c5v").textContent = fmt(fT) + " Hz";
	const t = [],
		v = [];
	for (let i = 0; i < N0; i++) {
		const ti = (i + jit[i]) / FS0;
		if (gap && ti >= 6 && ti < 6.3) continue;
		t.push(ti);
		v.push(
			aD * Math.sin(2 * Math.PI * 0.5 * ti + 1) +
				aT * Math.sin(2 * Math.PI * fT * ti) +
				aN * nz[i],
		);
	}
	const n = t.length;
	const fsEff = (n - 1) / (t[n - 1] - t[0]);
	let maxGap = 0;
	for (let i = 1; i < n; i++) maxGap = Math.max(maxGap, t[i] - t[i - 1]);
	const r = resample(t, v, 100);
	const xl = filtfilt(r.x, biquad("lp", 20, 100, 0.7071));
	let x50 = [];
	for (let i = 0; i < xl.length; i += 2) x50.push(xl[i]);
	const mean = x50.reduce((a, b) => a + b, 0) / x50.length;
	x50 = x50.map((a) => a - mean);
	const y = filtfilt(
		filtfilt(x50, biquad("hp", 1, 50, 0.7071)),
		biquad("lp", 20, 50, 0.7071),
	);
	const sy = y.slice(150),
		sx = x50.slice(150);
	const rms = (a) => Math.sqrt(a.reduce((s, q) => s + q * q, 0) / a.length);
	const rmsPre = rms(sx),
		rmsBand = rms(sy);
	const w = welch(sy);
	let fpk = 0,
		apk = 0,
		eb = 0,
		et = 0;
	const med = [];
	for (let k = 1; k <= 60; k++) {
		const f = k * 0.25,
			a = w.A[k];
		if (f >= 3) {
			med.push(a);
			if (a > apk) {
				apk = a;
				fpk = f;
			}
		}
		if (f >= 1) et += a * a;
		if (f >= 3.5 && f <= 7) eb += a * a;
	}
	med.sort((a, b) => a - b);
	const medv = med[Math.floor(med.length / 2)];
	const bpr = et > 0 ? eb / et : 0,
		snr = apk / Math.max(medv, 1e-6);
	let hit = 0;
	w.segAmp.forEach((a) => {
		let pk = 0,
			pf = 0;
		for (let k = 1; k <= 60; k++) {
			const f = k * 0.25;
			if (f >= 3 && a[k] > pk) {
				pk = a[k];
				pf = f;
			}
		}
		if (pf >= 3.5 && pf <= 7) hit++;
	});
	const cons = w.ns ? hit / w.ns : 0;
	const gates = [
		{ n: "Sampling di atas 50 Hz", ok: fsEff >= 50, d: fmt(fsEff, 0) + " Hz" },
		{
			n: "Tidak ada jeda lebih dari 200 ms",
			ok: maxGap <= 0.2,
			d: fmt(maxGap * 1000, 0) + " ms",
		},
		{
			n: "Gerakan tangan tidak terlalu besar",
			ok: rmsPre < 1.0,
			d: "RMS " + fmt(rmsPre, 2),
		},
		{
			n: "Sinyal di atas noise sensor",
			ok: rmsBand > 0.01,
			d: "RMS " + fmt(rmsBand, 3),
		},
	];
	const valid = gates.every((g) => g.ok);
	const detected =
		fpk >= 3.5 && fpk <= 7 && snr >= 4 && bpr >= 0.3 && cons >= 0.6;
	let ym = 0.5;
	for (let i = 0; i < x50.length; i++) ym = Math.max(ym, Math.abs(x50[i]));
	ym = Math.ceil(ym * 10) / 10;
	const t50 = x50.map((_, i) => r.t0 + i / 50);
	return {
		t: t,
		v: v,
		x50: x50,
		y: y,
		t50: t50,
		w: w,
		fpk: fpk,
		apk: apk,
		bpr: bpr,
		snr: snr,
		cons: cons,
		gates: gates,
		valid: valid,
		detected: detected,
		ym: ym,
	};
}

let stage = 1,
	S = null;
const X0 = 44,
	XW = 620,
	YC = 100,
	YH = 80;
const tx = (x, y, s, a) =>
	`<text x="${x}" y="${y}" text-anchor="${a || "middle"}" font-size="11" style="fill:var(--text-muted)">${s}</text>`;
function frame(w0, w1, ym) {
	let s = `<line x1="${X0}" y1="${YC}" x2="${X0 + XW}" y2="${YC}" style="stroke:var(--border-strong)" stroke-width="0.5"/>`;
	[w0, (w0 + w1) / 2, w1].forEach((q) => {
		const x = X0 + ((q - w0) / (w1 - w0)) * XW;
		s +=
			`<line x1="${x}" y1="${YC - YH}" x2="${x}" y2="${YC + YH}" style="stroke:var(--border)" stroke-width="0.5"/>` +
			tx(x, 200, fmt(q, 1) + " s");
	});
	s +=
		tx(X0 - 6, YC - YH + 8, "+" + fmt(ym, 1), "end") +
		tx(X0 - 6, YC + YH, "-" + fmt(ym, 1), "end");
	return s;
}
function pts(ts, vs, w0, w1, ym) {
	const o = [];
	for (let i = 0; i < ts.length; i++) {
		if (ts[i] < w0 || ts[i] > w1) continue;
		o.push([
			X0 + ((ts[i] - w0) / (w1 - w0)) * XW,
			YC - Math.max(-1, Math.min(1, vs[i] / ym)) * YH,
		]);
	}
	return o;
}
const path = (p) =>
	p
		.map((q, i) => (i ? "L" : "M") + q[0].toFixed(1) + " " + q[1].toFixed(1))
		.join("");
const dots = (p) =>
	p
		.map(
			(q) =>
				`<circle cx="${q[0].toFixed(1)}" cy="${q[1].toFixed(1)}" r="2.2" style="fill:var(--text-accent)"/>`,
		)
		.join("");

const caps = {
	1: "Titik asli dari sensor, 1 detik. Perhatikan jaraknya: ada yang rapat, ada yang renggang.",
	2: "Setelah resample: titik berjarak sama, 50 per detik.",
	3: "4 detik. Garis tipis abu-abu: sebelum filter. Garis tebal berwarna: sesudah filter, gerakan lambat hilang.",
	4: "Spektrum: tinggi batang adalah kekuatan tiap frekuensi. Area berwarna adalah pita tremor 3,5–7 Hz.",
};

function drawPlot() {
	let s = "";
	if (stage === 1) {
		const p = pts(S.t, S.v, 3, 4, S.ym);
		s =
			frame(3, 4, S.ym) +
			`<path d="${path(p)}" fill="none" style="stroke:var(--border-strong)" stroke-width="0.5"/>` +
			dots(p);
	} else if (stage === 2) {
		const p = pts(S.t50, S.x50, 3, 4, S.ym);
		s =
			frame(3, 4, S.ym) +
			`<path d="${path(p)}" fill="none" style="stroke:var(--border-strong)" stroke-width="0.5"/>` +
			dots(p);
	} else if (stage === 3) {
		const a = pts(S.t50, S.x50, 3, 7, S.ym),
			b = pts(S.t50, S.y, 3, 7, S.ym);
		s =
			frame(3, 7, S.ym) +
			`<path d="${path(a)}" fill="none" style="stroke:var(--text-muted)" stroke-width="1"/>` +
			`<path d="${path(b)}" fill="none" style="stroke:var(--text-accent)" stroke-width="1.8"/>`;
	} else {
		const base = 185,
			hh = 150,
			xf = (f) => X0 + (f / 15) * 600;
		let ym = 0.3;
		for (let k = 1; k <= 60; k++) ym = Math.max(ym, S.w.A[k]);
		s = `<rect x="${xf(3.5) - 3}" y="10" width="${xf(7) - xf(3.5) + 6}" height="${base - 10}" style="fill:var(--bg-accent)"/>`;
		s += tx((xf(3.5) + xf(7)) / 2, 24, "pita tremor");
		for (let k = 1; k <= 60; k++) {
			const f = k * 0.25,
				h = Math.min(1, S.w.A[k] / ym) * hh,
				inb = f >= 3.5 && f <= 7;
			s += `<rect x="${(xf(f) - 2.5).toFixed(1)}" y="${(base - h).toFixed(1)}" width="5" height="${h.toFixed(1)}" style="fill:${inb ? "var(--text-accent)" : "var(--text-muted)"}"/>`;
		}
		s += `<line x1="${X0}" y1="${base}" x2="${X0 + 600}" y2="${base}" style="stroke:var(--border-strong)" stroke-width="0.5"/>`;
		[0, 3, 6, 9, 12, 15].forEach((f) => {
			s += tx(xf(f), 200, f + " Hz");
		});
		if (S.apk > 0.02) {
			const h = Math.min(1, S.apk / ym) * hh;
			s += tx(xf(S.fpk), base - h - 6, fmt(S.fpk) + " Hz");
		}
	}
	$("plot").innerHTML = s;
	$("cap").textContent = caps[stage];
	document
		.querySelectorAll(".sb")
		.forEach((b) => b.classList.toggle("on", +b.dataset.s === stage));
}

function drawInfo() {
	$("m1").textContent = S.apk > 0.02 ? fmt(S.fpk) + " Hz" : "-";
	$("m2").textContent = fmt(S.bpr * 100, 0) + "%";
	$("m3").textContent = fmt(S.snr, 1) + " kali";
	$("m4").textContent = fmt(S.cons * 100, 0) + "%";
	$("gates").innerHTML = S.gates
		.map(
			(g) =>
				`<div class="g"><span>${g.n}</span><span style="color:var(${g.ok ? "--text-success" : "--text-danger"})">${g.ok ? "lolos" : "gagal"} (${g.d})</span></div>`,
		)
		.join("");
	const d = $("dec");
	if (!S.valid) {
		const why = S.gates
			.filter((g) => !g.ok)
			.map((g) => g.n.toLowerCase())
			.join("; ");
		d.style.background = "var(--bg-danger)";
		d.style.borderColor = "var(--border-danger)";
		d.style.color = "var(--text-danger)";
		d.textContent = "Pengukuran tidak valid, ulangi. Gagal pada: " + why + ".";
	} else if (S.detected) {
		d.style.background = "var(--bg-warning)";
		d.style.borderColor = "var(--border-warning)";
		d.style.color = "var(--text-warning)";
		d.textContent =
			"Tremor ritmis terdeteksi di pita khas Parkinson. Perlu evaluasi lebih lanjut.";
	} else {
		d.style.background = "var(--surface-1)";
		d.style.borderColor = "var(--border)";
		d.style.color = "var(--text-primary)";
		d.textContent = "Tremor ritmis tidak terdeteksi pada pengukuran ini.";
	}
}

function update() {
	S = proc();
	drawPlot();
	drawInfo();
}
["c1", "c2", "c3", "c5"].forEach((i) => $(i).addEventListener("input", update));
$("c4").addEventListener("change", update);
document.querySelectorAll(".sb").forEach((b) =>
	b.addEventListener("click", () => {
		stage = +b.dataset.s;
		drawPlot();
	}),
);
update();
