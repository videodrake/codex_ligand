import path from "node:path";
import fs from "node:fs/promises";
import { fileURLToPath } from "node:url";
import {
  Presentation,
  PresentationFile,
  layers,
  column,
  row,
  grid,
  panel,
  text,
  image,
  shape,
  rule,
  fill,
  hug,
  fixed,
  wrap,
  fr,
} from "@oai/artifact-tool";

const root = "C:/Users/dudwn/OneDrive/문서/EGFR_MYO1D/analysis_v0_5";
const assets = path.join(root, "reports/assets");
const outDir = path.join(root, "reports/ppt_output");
const previewDir = path.join(outDir, "previews");
const untilArg = process.argv.find((arg) => arg.startsWith("--until="));
const until = untilArg ? Number(untilArg.split("=")[1]) : 1;

const W = 1920;
const H = 1080;
const C = {
  ink: "#132026",
  slate: "#40505A",
  muted: "#65747C",
  pale: "#EFF7F8",
  cyan: "#00A8C8",
  green: "#23A66D",
  orange: "#E6903A",
  navy: "#0B2630",
  white: "#FFFFFF",
  line: "#D7E4E8",
  soft: "#F5FAFB",
};

function asset(name) {
  return path.join(assets, name).replaceAll("\\", "/");
}

async function imageDataUrl(filePath) {
  const ext = path.extname(filePath).toLowerCase();
  const mime = ext === ".jpg" || ext === ".jpeg" ? "image/jpeg" : "image/png";
  const bytes = await fs.readFile(filePath);
  return `data:${mime};base64,${bytes.toString("base64")}`;
}

async function optionalImageDataUrl(filePath) {
  try {
    return await imageDataUrl(filePath);
  } catch {
    return null;
  }
}

const coverBackgroundDataUrl = await imageDataUrl(asset("slide01_cover_background.png"));
const slide04RmsdDataUrl = await imageDataUrl(asset("slide04_rmsd_by_orientation.png"));
const slide04RgDataUrl = await imageDataUrl(asset("slide04_rg_by_orientation.png"));
const slide04HbondDataUrl = await imageDataUrl(asset("slide04_hbond_by_orientation.png"));
const slide04ScoreDataUrl = await imageDataUrl(asset("slide04_orientation_score_summary.png"));
const slide05KdDistanceDataUrl = await imageDataUrl(asset("slide05_kd_membrane_distance_200ns.png"));
const slide05ContactsDataUrl = await imageDataUrl(asset("slide05_interchain_contacts_200ns.png"));
const slide05TiltDataUrl = await imageDataUrl(asset("slide05_tm_kd_tilt_200ns.png"));
const slide06SideDataUrl = await optionalImageDataUrl(asset("slide6_plus10_side_view.png"));
const slide06TopDataUrl = await optionalImageDataUrl(asset("slide6_plus10_top_view.png"));
const selectedClusterDir = path.join(root, "3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173");
const slide07ClusterPopulationDataUrl = await imageDataUrl(path.join(selectedClusterDir, "cluster_population.png"));
const slide07ClusterTimelineDataUrl = await imageDataUrl(path.join(selectedClusterDir, "cluster_timeline.png"));

function addFooter(slide, label, source = "") {
  slide.compose(
    row({ name: `footer-${label}`, width: fill, height: hug, padding: { x: 84, y: 32 }, align: "center", justify: "between" }, [
      text(label, { name: `footer-label-${label}`, width: fixed(90), height: hug, style: { fontSize: 15, color: C.muted } }),
      text(source, { name: `footer-source-${label}`, width: wrap(1100), height: hug, style: { fontSize: 12, color: "#7C8A91" } }),
    ]),
    { frame: { left: 0, top: 970, width: W, height: 100 }, baseUnit: 8 },
  );
}

function eyebrow(value) {
  return text(value, { name: "eyebrow", width: fixed(560), height: hug, style: { fontSize: 20, bold: true, color: C.cyan } });
}

function slideTitle(value, opts = {}) {
  return text(value, {
    name: opts.name ?? "slide-title",
    width: opts.width ?? wrap(1420),
    height: hug,
    style: { fontSize: opts.fontSize ?? 52, bold: true, color: opts.color ?? C.ink },
  });
}

function subtitle(value, opts = {}) {
  return text(value, {
    name: opts.name ?? "slide-subtitle",
    width: opts.width ?? wrap(1220),
    height: hug,
    style: { fontSize: opts.fontSize ?? 26, color: opts.color ?? C.slate },
  });
}

function smallCaps(value, color = C.cyan) {
  return text(value, { width: wrap(520), height: hug, style: { fontSize: 16, bold: true, color } });
}

function pillLabel(value, color = C.cyan) {
  return text(value, { width: fill, height: hug, style: { fontSize: 14, bold: true, color, alignment: "center" } });
}

function bullet(value, color = C.ink) {
  return row({ width: fill, height: hug, gap: 14, align: "start" }, [
    shape({ width: fixed(9), height: fixed(9), fill: C.cyan, borderRadius: "rounded-full" }),
    text(value, { width: fill, height: hug, style: { fontSize: 23, color } }),
  ]);
}

function metric(label, value, note, accent = C.cyan) {
  return panel({ width: fill, height: fixed(155), padding: { x: 24, y: 20 }, fill: C.white, borderRadius: 10 },
    column({ width: fill, height: fill, gap: 6 }, [
      text(label, { width: fill, height: hug, style: { fontSize: 15, bold: true, color: C.muted } }),
      text(value, { width: fill, height: hug, style: { fontSize: 42, bold: true, color: accent } }),
      text(note, { width: fill, height: hug, style: { fontSize: 14, color: C.muted } }),
    ]),
  );
}

function stepNode(index, title, note, status, accent = C.cyan) {
  const fillColor = status === "current" ? "#E8F7F2" : status === "future" ? "#F8FBFC" : "#E8F5FA";
  return panel({ width: fixed(270), height: fixed(174), padding: { x: 22, y: 18 }, fill: fillColor, borderRadius: 8 },
    column({ width: fill, height: fill, gap: 10 }, [
      text(String(index).padStart(2, "0"), { width: fill, height: hug, style: { fontSize: 18, bold: true, color: accent } }),
      text(title, { width: fill, height: hug, style: { fontSize: 23, bold: true, color: C.ink } }),
      text(note, { width: fill, height: hug, style: { fontSize: 16, color: C.muted } }),
    ]),
  );
}

function arrowMark() {
  return text("→", { width: fixed(42), height: fixed(174), style: { fontSize: 34, bold: true, color: "#9AB1BA", alignment: "center" } });
}

function tableCell(value, width, opts = {}) {
  return text(value, {
    width: fixed(width),
    height: hug,
    style: { fontSize: opts.fontSize ?? 15, bold: opts.bold ?? false, color: opts.color ?? C.ink },
  });
}

function templateRow(component, source, role, caveat, fillColor = "#F8FBFC") {
  return panel({ width: fill, height: fixed(82), padding: { x: 18, y: 13 }, fill: fillColor, borderRadius: 6 },
    row({ width: fill, height: fill, gap: 18, align: "center" }, [
      tableCell(component, 180, { bold: true, fontSize: 16 }),
      tableCell(source, 130, { bold: true, color: C.cyan, fontSize: 16 }),
      tableCell(role, 355, { color: C.slate }),
      tableCell(caveat, 330, { color: C.muted, fontSize: 14 }),
    ]),
  );
}

function workflowStep(title, note, accent = C.cyan) {
  return panel({ width: fill, height: fixed(74), padding: { x: 20, y: 13 }, fill: "#F8FBFC", borderRadius: 8 },
    row({ width: fill, height: fill, gap: 16, align: "center" }, [
      shape({ width: fixed(8), height: fixed(42), fill: accent, borderRadius: "rounded-full" }),
      column({ width: fill, height: hug, gap: 4 }, [
        text(title, { width: fill, height: hug, style: { fontSize: 19, bold: true, color: C.ink } }),
        text(note, { width: fill, height: hug, style: { fontSize: 14, color: C.muted } }),
      ]),
    ]),
  );
}

function figurePanel(name, dataUrl, caption) {
  return panel({ width: fill, height: fill, padding: { x: 12, y: 10 }, fill: "#F8FBFC", borderRadius: 8 },
    column({ width: fill, height: fill, gap: 6 }, [
      image({ name, dataUrl, contentType: "image/png", width: fill, height: fill, fit: "contain", alt: caption }),
      text(caption, { width: fill, height: hug, style: { fontSize: 12, color: C.muted } }),
    ]),
  );
}

function renderSlot(name, dataUrl, title, outputPath, sourcePath, scriptPath) {
  if (dataUrl) {
    return figurePanel(name, dataUrl, title);
  }
  return panel({ width: fill, height: fill, padding: { x: 28, y: 24 }, fill: "#F8FBFC", borderRadius: 8 },
    column({ width: fill, height: fill, gap: 14 }, [
      smallCaps("render image needed", C.orange),
      text(title, { width: fill, height: hug, style: { fontSize: 28, bold: true, color: C.ink } }),
      rule({ width: fixed(160), stroke: C.orange, weight: 4 }),
      text(`Add later: ${outputPath}`, { width: fill, height: hug, style: { fontSize: 16, color: C.slate } }),
      text(`Source: ${sourcePath}`, { width: fill, height: hug, style: { fontSize: 14, color: C.muted } }),
      text(`Render script: ${scriptPath}`, { width: fill, height: hug, style: { fontSize: 14, color: C.muted } }),
    ]),
  );
}

function statusRow(item, finding, status, accent = C.orange) {
  return panel({ width: fill, height: fixed(78), padding: { x: 18, y: 12 }, fill: "#F8FBFC", borderRadius: 6 },
    row({ width: fill, height: fill, gap: 14, align: "center" }, [
      tableCell(item, 255, { bold: true, fontSize: 15 }),
      tableCell(finding, 330, { color: C.slate, fontSize: 14 }),
      text(status, { width: fixed(150), height: hug, style: { fontSize: 14, bold: true, color: accent } }),
    ]),
  );
}

function statusCard(title, items, accent = C.cyan, fillColor = "#F8FBFC") {
  return panel({ width: fill, height: fill, padding: { x: 24, y: 22 }, fill: fillColor, borderRadius: 8 },
    column({ width: fill, height: fill, gap: 18 }, [
      smallCaps(title, accent),
      ...items.map((item) => bullet(item, C.ink)),
    ]),
  );
}

function roadmapStep(index, title, note, accent = C.cyan) {
  return row({ width: fill, height: fixed(74), gap: 16, align: "center" }, [
    panel({ width: fixed(58), height: fixed(58), padding: { x: 14, y: 14 }, fill: accent, borderRadius: "rounded-full" },
      text(String(index), { width: fill, height: hug, style: { fontSize: 20, bold: true, color: C.white, alignment: "center" } }),
    ),
    column({ width: fill, height: hug, gap: 4 }, [
      text(title, { width: fill, height: hug, style: { fontSize: 22, bold: true, color: C.ink } }),
      text(note, { width: fill, height: hug, style: { fontSize: 15, color: C.muted } }),
    ]),
  ]);
}

const presentation = Presentation.create({ slideSize: { width: W, height: H } });

function makeSlide1() {
  const slide = presentation.slides.add();
  slide.compose(
    image({ name: "cover-background", dataUrl: coverBackgroundDataUrl, contentType: "image/png", width: fill, height: fill, fit: "cover", alt: "Abstract scientific membrane receptor background" }),
    { frame: { left: 0, top: 0, width: W, height: H }, baseUnit: 8 },
  );
  slide.compose(
    column({ name: "cover-text-stack", width: fill, height: hug, gap: 19 }, [
      eyebrow("EGFR-MYO1D 진행 보고"),
      column({ name: "title-lines", width: fill, height: hug, gap: 8 }, [
        text("EGFR 모델을 먼저", { name: "slide-title-line-1", width: fill, height: hug, style: { fontSize: 62, bold: true, color: C.ink } }),
        text("만들고 확인했습니다", { name: "slide-title-line-2", width: fill, height: hug, style: { fontSize: 62, bold: true, color: C.ink } }),
        text("대표 구조도 정했습니다", { name: "slide-title-line-3", width: fill, height: hug, style: { fontSize: 58, bold: true, color: C.ink } }),
      ]),
      rule({ width: fixed(250), stroke: C.green, weight: 5 }),
      subtitle("MYO1D가 EGFR에 붙는지 보려면, 먼저 EGFR의 모양과 방향을 정해야 합니다.", { width: wrap(900), fontSize: 25 }),
      row({ width: fixed(930), height: hug, gap: 14, padding: { y: 8 }, align: "center" }, [
        panel({ width: fixed(280), height: fixed(44), padding: { x: 18, y: 10 }, fill: "#E8F7F2", borderRadius: "rounded-full" }, pillLabel("모델 만들기", C.green)),
        panel({ width: fixed(230), height: fixed(44), padding: { x: 18, y: 10 }, fill: "#E8F5FA", borderRadius: "rounded-full" }, pillLabel("200 ns 확인", C.cyan)),
        panel({ width: fixed(220), height: fixed(44), padding: { x: 18, y: 10 }, fill: "#FFF2E4", borderRadius: "rounded-full" }, pillLabel("대표 구조 선택", C.orange)),
      ]),
      text("이번 발표는 EGFR 대표 구조를 정하고 MYO1D 분석으로 넘어가는 준비 상황입니다.", { name: "scope-note", width: wrap(930), height: hug, style: { fontSize: 17, color: C.muted } }),
    ]),
    { frame: { left: 86, top: 76, width: 980, height: 790 }, baseUnit: 8 },
  );
  addFooter(slide, "01 / 09", "표지 그림은 개념 그림입니다. 결과 설명은 로컬 데이터 파일을 기준으로 했습니다.");
}

function makeSlide2() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide2-root", width: fill, height: fill, padding: { x: 84, y: 68 }, gap: 30 }, [
      column({ name: "slide2-title", width: fill, height: hug, gap: 12 }, [
        eyebrow("왜 먼저 해야 하나?"),
        slideTitle("MYO1D를 보기 전에 EGFR 모양부터 정합니다"),
        subtitle("EGFR 모양이 바뀌면 MYO1D가 붙는 위치도 달라 보일 수 있습니다. 그래서 기준이 될 EGFR 모델을 먼저 정했습니다.", { width: wrap(1320) }),
      ]),
      row({ name: "slide2-flow", width: fill, height: fixed(188), gap: 18, align: "center" }, [
        stepNode(1, "EGFR 대표 구조 정하기", "Cluster 1/2 구간 대표 구조 사용", "current", C.green),
        arrowMark(),
        stepNode(2, "MYO1D 붙여보기", "EGFR 모델이 정해진 뒤 진행", "future", C.cyan),
        arrowMark(),
        stepNode(3, "붙는 자리 찾기", "EGFR의 어느 부분과 닿는지 확인", "future", C.cyan),
        arrowMark(),
        stepNode(4, "작은 홈 찾기", "약물이 들어갈 만한 공간 후보", "future", C.orange),
        arrowMark(),
        stepNode(5, "화합물 넣어보기", "이 발표의 결과가 아니라 다음 단계", "future", C.orange),
      ]),
      row({ name: "slide2-bottom", width: fill, height: fill, gap: 48 }, [
        column({ name: "slide2-confirmed", width: fill, height: fill, gap: 16 }, [
          smallCaps("이번 발표 내용", C.green),
          bullet("EGFR의 여러 부분을 이어 붙여 기본 모델을 만들었습니다."),
          bullet("+10° 방향 모델이 현재 지표에서 가장 좋아 보였습니다."),
          bullet("Cluster 1/2 대표 구조를 다음 입력 후보로 정했습니다."),
        ]),
        column({ name: "slide2-boundary", width: fill, height: fill, gap: 16 }, [
          smallCaps("아직 다음 단계", C.orange),
          bullet("MYO1D를 붙여본 최종 결과는 아직 아닙니다."),
          bullet("붙는 자리와 작은 홈 찾기는 다음 단계입니다."),
          bullet("화합물 docking은 MYO1D 자리와 pocket을 본 뒤 시작합니다."),
        ]),
      ]),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "02 / 09", "자료: reports/EGFR_modeling_clustering_PPT_data_package.md");
}

function makeSlide3() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide3-root", width: fill, height: fill, padding: { x: 84, y: 64 }, gap: 24 }, [
      column({ name: "slide3-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("모델을 만든 방법"),
        slideTitle("EGFR은 여러 자료를 이어 붙여 만들었습니다"),
        subtitle("EGFR의 막 부분, 연결 부분, 세포 안쪽 효소 부분을 각각 다른 구조 자료에서 가져와 하나의 모델로 만들었습니다.", { width: wrap(1480), fontSize: 24 }),
      ]),
      row({ name: "slide3-body", width: fill, height: fill, gap: 44 }, [
        column({ name: "slide3-template-table", width: fixed(1110), height: fill, gap: 7 }, [
          panel({ width: fill, height: fixed(50), padding: { x: 18, y: 12 }, fill: C.navy, borderRadius: 6 },
            row({ width: fill, height: fill, gap: 18, align: "center" }, [
              tableCell("부분", 180, { bold: true, color: C.white, fontSize: 15 }),
              tableCell("자료", 130, { bold: true, color: C.white, fontSize: 15 }),
              tableCell("무엇에 썼나", 355, { bold: true, color: C.white, fontSize: 15 }),
              tableCell("주의할 점", 330, { bold: true, color: C.white, fontSize: 15 }),
            ]),
          ),
          templateRow("막을 지나는 부분", "2M0B", "EGFR 두 가닥이 막을 지나는 모양", "스크립트에서는 634-670 구간을 사용", "#E8F7F2"),
          templateRow("막 근처 연결부", "2M20", "막 근처에서 이어지는 방향 참고", "좌표 생성 과정 일부는 추가 확인 필요"),
          templateRow("세포 안쪽 효소부", "3GT8", "inactive KD dimer 모양 참고", "번호 체계가 달라서 해석 시 주의"),
          templateRow("빈 구간 채우기", "MODELLER", "없는 구간을 계산으로 메움", "실험 구조가 아니라 모델링 결과", "#FFF8EE"),
        ]),
        column({ name: "slide3-workflow", width: fill, height: fill, gap: 12 }, [
          smallCaps("작업 순서", C.green),
          workflowStep("막 부분을 먼저 맞춤", "2M0B 자료로 EGFR 두 가닥의 막 부분을 잡음", C.green),
          workflowStep("연결부 방향을 참고", "2M20 자료로 막 근처 방향을 정리", C.cyan),
          workflowStep("세포 안쪽 큰 부분을 붙임", "3GT8 자료로 KD dimer 부분을 배치", C.cyan),
          workflowStep("빈 구간을 계산으로 채움", "JM-B와 A-loop처럼 없는 부분을 MODELLER로 보완", C.orange),
          workflowStep("방향 후보 5개를 만듦", "-20°, -10°, 0°, +10°, +20°를 비교", C.green),
          text("Slide 6의 구조 그림은 이 모델 PDB의 실제 CA 좌표를 이용해 만든 overview입니다.", { width: fill, height: hug, style: { fontSize: 14, color: C.muted } }),
        ]),
      ]),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "03 / 09", "자료: 1.align/*, 01_numbering_provenance/template_model_mapping.csv");
}

function makeSlide4() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide4-root", width: fill, height: fill, padding: { x: 84, y: 62 }, gap: 24 }, [
      column({ name: "slide4-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("방향 비교"),
        slideTitle("+10° 방향 모델이 가장 좋아 보였습니다"),
        subtitle("EGFR 안쪽 큰 부분을 여러 각도로 돌려 비교했고, 현재 계산된 지표에서는 +10°가 가장 안정적으로 보였습니다.", { width: wrap(1460), fontSize: 24 }),
      ]),
      row({ name: "slide4-body", width: fill, height: fill, gap: 40 }, [
        column({ name: "slide4-metrics", width: fixed(500), height: fill, gap: 10 }, [
          smallCaps("선택 이유", C.green),
          metric("움직임 크기", "2.68 Å", "RMSD가 낮을수록 덜 흔들린다는 뜻", C.green),
          metric("전체 크기", "28.31 Å", "Rg가 비슷하게 유지되어 비교적 안정적", C.cyan),
          metric("두 가닥 연결", "11.49", "H-bond가 더 많이 보인 조건", C.orange),
          text("주의: 모든 지표가 완벽하게 다 있는 것은 아닙니다. 그래도 현재 가지고 있는 자료에서는 +10°가 가장 좋은 후보입니다.", { width: fill, height: hug, style: { fontSize: 14, color: C.muted } }),
        ]),
        grid({ name: "slide4-figures", width: fill, height: fill, columns: [fr(1), fr(1)], rows: [fr(1), fr(1)], columnGap: 18, rowGap: 18 }, [
          figurePanel("slide4-rmsd", slide04RmsdDataUrl, "각 방향별 흔들림 정도(RMSD)"),
          figurePanel("slide4-rg", slide04RgDataUrl, "각 방향별 전체 크기 변화(Rg)"),
          figurePanel("slide4-hbond", slide04HbondDataUrl, "두 가닥 사이 연결 정도(H-bond)"),
          figurePanel("slide4-score", slide04ScoreDataUrl, "현재 지표를 합친 방향 점수"),
        ]),
      ]),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "04 / 09", "자료: 02_orientation_metrics/orientation_comparison_table.csv, figures/");
}

function makeSlide5() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide5-root", width: fill, height: fill, padding: { x: 84, y: 60 }, gap: 22 }, [
      column({ name: "slide5-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("200 ns 확인"),
        slideTitle("최종 움직임 기록 파일을 확인했습니다"),
        subtitle("trajectory는 단백질이 시간에 따라 어떻게 움직였는지 저장한 기록입니다. 이번에는 +10° 모델의 200 ns 기록을 확인했습니다.", { width: wrap(1480), fontSize: 24 }),
      ]),
      row({ name: "slide5-gmx-strip", width: fill, height: fixed(155), gap: 18 }, [
        metric("원자 수", "49,854", "물 분자를 뺀 최종 시스템", C.green),
        metric("저장 장면 수", "2,001", "시간별로 저장된 구조 장면", C.cyan),
        metric("확인한 시간", "0-200 ns", "ns는 아주 짧은 시간 단위", C.orange),
        metric("검사 결과", "OK", "분석에 사용할 수 있는 상태", C.green),
      ]),
      grid({ name: "slide5-figures", width: fill, height: fill, columns: [fr(1), fr(1), fr(1)], rows: [fr(1)], columnGap: 18 }, [
        figurePanel("slide5-kd-distance", slide05KdDistanceDataUrl, "KD 부분이 막에서 얼마나 떨어져 있는지"),
        figurePanel("slide5-contacts", slide05ContactsDataUrl, "EGFR 두 가닥이 서로 얼마나 닿아 있는지"),
        figurePanel("slide5-tm-kd-tilt", slide05TiltDataUrl, "막 부분과 KD 부분의 각도 변화"),
      ]),
      text("사용한 최종 파일: EGFR_plus10_step7_production_nw.gro + EGFR_plus10_step7_200ns_nw.xtc. 다른 full_200ns.xtc 파일은 최종 근거로 쓰지 않았습니다.", { width: fill, height: hug, style: { fontSize: 15, color: C.muted } }),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "05 / 09", "자료: gmx_check_EGFR_plus10_step7_200ns_nw.txt, trajectory_metrics/");
}

function makeSlide6() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide6-root", width: fill, height: fill, padding: { x: 84, y: 62 }, gap: 22 }, [
      column({ name: "slide6-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("구조 그림으로 확인"),
        slideTitle("EGFR 구조 그림", { width: fixed(900), fontSize: 46 }),
        subtitle("아래 그림은 실제 PDB 좌표를 간단히 그린 overview입니다. 현재 진행 상황 설명용으로 사용합니다.", { width: wrap(1480), fontSize: 24 }),
      ]),
      row({ name: "slide6-body", width: fill, height: fill, gap: 34 }, [
        grid({ name: "slide6-render-slots", width: fixed(1080), height: fill, columns: [fr(1), fr(1)], rows: [fr(1)], columnGap: 18 }, [
          renderSlot("slide6-side", slide06SideDataUrl, "옆에서 본 EGFR 모델과 막 위치", "reports/assets/slide6_plus10_side_view.png", "03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb", "reports/ppt_src/make_structure_overview.py"),
          renderSlot("slide6-top", slide06TopDataUrl, "위에서 본 EGFR 두 가닥의 배치", "reports/assets/slide6_plus10_top_view.png", "03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb", "reports/ppt_src/make_structure_overview.py"),
        ]),
        column({ name: "slide6-notes", width: fill, height: fill, gap: 16 }, [
          smallCaps("쉽게 말하면", C.green),
          bullet("+10° 모델은 막 근처에서 EGFR 두 가닥이 같이 있는 모양을 유지합니다."),
          bullet("200 ns 확인 결과, 다음 MYO1D 분석에 넣을 후보로 쓸 수 있어 보입니다."),
          metric("막과의 거리", "25.66 Å", "150-200 ns 구간 평균", C.green),
          metric("두 가닥 접촉", "505.28", "서로 가까운 원자 접촉 수", C.cyan),
          text("주의: 이 그림은 실제 좌표를 단순하게 그린 overview입니다. 더 예쁜 최종 그림은 PyMOL/ChimeraX로 나중에 만들 수 있습니다.", { width: fill, height: hug, style: { fontSize: 15, color: C.muted } }),
        ]),
      ]),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "06 / 09", "자료: EGFR_plus10_final_model_EGFR_rot10_best.pdb, membrane_stability_*_mdanalysis.csv");
}

function makeSlide7() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide7-root", width: fill, height: fill, padding: { x: 84, y: 62 }, gap: 22 }, [
      column({ name: "slide7-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("대표 구조 선택"),
        slideTitle("Cluster 1과 2 대표 구조를 쓰기로 했습니다"),
        subtitle("clustering은 비슷한 구조끼리 묶어서 대표 모양을 고르는 과정입니다. 이번에는 Cluster 1/2가 나오는 구간을 다음 MYO1D 분석 입력으로 사용합니다.", { width: wrap(1510), fontSize: 24 }),
      ]),
      row({ name: "slide7-body", width: fill, height: fill, gap: 36 }, [
        grid({ name: "slide7-figures", width: fixed(850), height: fill, columns: [fr(1)], rows: [fr(1), fr(1)], rowGap: 18 }, [
          figurePanel("slide7-population", slide07ClusterPopulationDataUrl, "선택 구간에서 C1과 C2가 가장 큰 대표 그룹으로 나타남"),
          figurePanel("slide7-timeline", slide07ClusterTimelineDataUrl, "시간에 따라 어느 cluster에 속했는지 표시"),
        ]),
        column({ name: "slide7-status", width: fill, height: fill, gap: 10 }, [
          smallCaps("선택한 대표 구조", C.green),
          statusRow("선택 구간", "70-100 ns, cutoff 0.173 nm", "사용", C.green),
          statusRow("Cluster 1", "121/301장, 40.2%, 대표 시간 79.3 ns", "대표 1", C.green),
          statusRow("Cluster 2", "68/301장, 22.6%, 대표 시간 93.5 ns", "대표 2", C.cyan),
          statusRow("C1 구조 파일", "EGFR_170-200.pdb로 정리된 대표 구조", "사용", C.green),
          statusRow("C2 구조 파일", "cluster_2_representative.pdb", "사용", C.cyan),
          text("쉽게 말하면: EGFR 대표 구조는 C1 하나만 보지 않고, C1과 C2 두 모양을 다음 MYO1D docking 입력 후보로 가져갑니다.", { width: fill, height: hug, style: { fontSize: 16, color: C.muted } }),
        ]),
      ]),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "07 / 09", "자료: 3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/clustering_summary.txt");
}

function makeSlide8() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide8-root", width: fill, height: fill, padding: { x: 84, y: 62 }, gap: 28 }, [
      column({ name: "slide8-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("현재 상태 정리"),
        slideTitle("완료한 일과 남은 일"),
        subtitle("발표에서는 완료된 일, 확인 중인 일, 다음에 할 일을 나누어 보여주면 됩니다.", { width: wrap(1500), fontSize: 24 }),
      ]),
      row({ name: "slide8-lanes", width: fill, height: fill, gap: 22 }, [
        statusCard("지금 말할 수 있음", [
          "EGFR 기본 모델을 만들었습니다.",
          "+10° 방향 모델을 현재 후보로 골랐습니다.",
          "200 ns 움직임 기록 파일을 확인했습니다.",
          "Cluster 1/2 대표 구조를 다음 입력 후보로 정했습니다.",
        ], C.green, "#E8F7F2"),
        statusCard("아직 확인 필요", [
          "C1/C2 구조 파일을 MYO1D docking 입력 형식으로 정리합니다.",
          "대표 구조별로 MYO1D docking 결과를 비교해야 합니다.",
          "더 예쁜 구조 그림은 나중에 PyMOL/ChimeraX로 만들 수 있습니다.",
          "일부 연결부 제작 과정은 추가 확인이 필요합니다.",
        ], C.orange, "#FFF8EE"),
        statusCard("다음 단계", [
          "MYO1D를 EGFR에 붙여보는 docking.",
          "MYO1D가 닿는 EGFR 위치 찾기.",
          "약물이 들어갈 만한 작은 홈 찾기.",
          "대표 구조가 정해진 뒤 화합물 docking.",
        ], C.cyan, "#E8F5FA"),
      ]),
      panel({ width: fill, height: fixed(74), padding: { x: 24, y: 16 }, fill: C.navy, borderRadius: 8 },
        text("발표 핵심: EGFR 모델 준비, 200 ns 확인, C1/C2 대표 구조 선택까지 진행 완료. 이제 MYO1D 분석으로 넘어갑니다.", { width: fill, height: hug, style: { fontSize: 20, bold: true, color: C.white } }),
      ),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "08 / 09", "자료: reports/EGFR_modeling_clustering_PPT_data_package.md, reports/다음할일.md");
}

function makeSlide9() {
  const slide = presentation.slides.add();
  slide.compose(
    column({ name: "slide9-root", width: fill, height: fill, padding: { x: 84, y: 62 }, gap: 28 }, [
      column({ name: "slide9-title", width: fill, height: hug, gap: 10 }, [
        eyebrow("다음 단계"),
        slideTitle("대표 구조를 정했고 MYO1D로 갑니다"),
        subtitle("다음 분석에서는 Cluster 1과 Cluster 2 대표 구조에 MYO1D를 붙여보고, 어느 위치가 반복해서 나오는지 비교합니다.", { width: wrap(1480), fontSize: 24 }),
      ]),
      row({ name: "slide9-body", width: fill, height: fill, gap: 44 }, [
        panel({ width: fixed(770), height: fill, padding: { x: 34, y: 30 }, fill: "#F8FBFC", borderRadius: 8 },
          column({ width: fill, height: fill, gap: 22 }, [
            smallCaps("정해진 입력", C.green),
            text("C1과 C2 대표 구조를 MYO1D docking에 사용합니다", { width: fill, height: hug, style: { fontSize: 36, bold: true, color: C.ink } }),
            rule({ width: fixed(200), stroke: C.orange, weight: 5 }),
            bullet("Cluster 1은 가장 많이 나온 대표 모양입니다."),
            bullet("Cluster 2는 두 번째로 많이 나온 대표 모양입니다."),
            bullet("689/692/713/715 같은 residue 이름을 최종 topology 기준으로 다시 확인합니다."),
            text("다음에 추가할 자료: C1/C2별 MYO1D docking 결과, 반복해서 닿는 EGFR 위치, pocket 후보.", { width: fill, height: hug, style: { fontSize: 16, color: C.muted } }),
          ]),
        ),
        column({ name: "slide9-roadmap", width: fill, height: fill, gap: 20 }, [
          smallCaps("다음 작업 순서", C.green),
          roadmapStep(1, "C1/C2 구조 준비", "두 대표 구조를 docking 입력으로 정리", C.green),
          roadmapStep(2, "MYO1D 붙여보기", "C1과 C2 각각에 MYO1D docking", C.cyan),
          roadmapStep(3, "닿는 자리 정리하기", "EGFR에서 MYO1D가 닿는 부분 표시", C.cyan),
          roadmapStep(4, "작은 홈 찾기", "약물이 들어갈 만한 공간 후보 찾기", C.orange),
          roadmapStep(5, "화합물 넣어보기", "대표 구조와 pocket이 정해진 뒤 진행", C.orange),
        ]),
      ]),
    ]),
    { frame: { left: 0, top: 0, width: W, height: 970 }, baseUnit: 8 },
  );
  addFooter(slide, "09 / 09", "다음 작업은 reports/다음할일.md에 따로 정리했습니다.");
}

const slideBuilders = [makeSlide1, makeSlide2, makeSlide3, makeSlide4, makeSlide5, makeSlide6, makeSlide7, makeSlide8, makeSlide9];
const builtSlides = Math.min(until, slideBuilders.length);
for (let i = 0; i < builtSlides; i++) slideBuilders[i]();

await fs.mkdir(outDir, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

const pptxBlob = await PresentationFile.exportPptx(presentation);
const progressPptx = path.join(outDir, `EGFR_MYO1D_progress_until_slide${builtSlides}.pptx`);
const legacyFinalPptx = path.join(outDir, "EGFR_MYO1D_receptor_modeling_clustering_deck.pptx");
const finalPptx = path.join(outDir, "EGFR_MYO1D_receptor_modeling_clustering_deck_cluster12_selected.pptx");
await pptxBlob.save(progressPptx);
if (builtSlides === slideBuilders.length) {
  await pptxBlob.save(finalPptx);
  try {
    await pptxBlob.save(legacyFinalPptx);
  } catch (error) {
    console.warn(JSON.stringify({ legacyFinalPptx, saveStatus: "skipped", reason: error?.code ?? String(error) }));
  }
}

for (let i = 0; i < presentation.slides.items.length; i++) {
  const png = await presentation.slides.items[i].export({ format: "png" });
  await fs.writeFile(path.join(previewDir, `slide${String(i + 1).padStart(2, "0")}.png`), Buffer.from(await png.arrayBuffer()));
}
console.log(JSON.stringify({ slides: presentation.slides.items.length, pptx: progressPptx, finalPptx: builtSlides === slideBuilders.length ? finalPptx : null, previewDir }, null, 2));
