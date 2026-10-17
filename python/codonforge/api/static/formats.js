/* FASTA, GenBank, and CSV for a pareto_front document.

   The bytes match python/codonforge/cli/writers.py. This file renders the
   document. It does not recompute a score. */

export function formatNumber(value) {
  if (value === null || value === undefined) return "";
  return Number(value).toFixed(6);
}

export function toFasta(front) {
  const order = front.objective_order;
  const sense = front.objective_sense;
  const lines = [];
  for (const candidate of front.candidates) {
    const values = order
      .map(
        (term, index) =>
          `${term}=${formatNumber(candidate.objective_vector[index])}(${sense[index]})`,
      )
      .join(" ");
    const feasible = candidate.constraint_report.feasible ? "feasible" : "INFEASIBLE";
    lines.push(
      `>${candidate.candidate_id} job=${front.job_id} origin=${candidate.provenance.origin} ${feasible} ${values}`,
    );
    const sequence = candidate.sequences[0].cds;
    for (let start = 0; start < sequence.length; start += 60) {
      lines.push(sequence.slice(start, start + 60));
    }
  }
  return `${lines.join("\n")}\n`;
}

function csvCell(value) {
  if (value === true) return "True";
  if (value === false) return "False";
  if (value === null || value === undefined) return "";
  const text = String(value);
  if (/[",\n\r]/.test(text)) return `"${text.replaceAll('"', '""')}"`;
  return text;
}

export function toCsv(front) {
  const header = ["candidate_id", "origin", "feasible", "nondomination_rank"];
  for (let index = 0; index < front.objective_order.length; index += 1) {
    header.push(`${front.objective_order[index]} (${front.objective_sense[index]})`);
  }
  header.push("cds");
  const lines = [header.map(csvCell).join(",")];
  for (const candidate of front.candidates) {
    const row = [
      candidate.candidate_id,
      candidate.provenance.origin,
      candidate.constraint_report.feasible,
      candidate.nondomination_rank,
    ];
    for (const value of candidate.objective_vector) row.push(formatNumber(value));
    row.push(candidate.sequences[0].cds);
    lines.push(row.map(csvCell).join(","));
  }
  return `${lines.join("\n")}\n`;
}

export function toGenbank(front, sequenceDocument, candidateIndex = 0) {
  const candidate = front.candidates[candidateIndex];
  const cds = sequenceDocument.cds;
  const five = sequenceDocument.five_prime_flank || "";
  const three = sequenceDocument.three_prime_flank || "";
  const construct = `${five}${cds}${three}`;
  const locus = candidate.candidate_id.slice(0, 16).padEnd(16, " ");
  const lines = [
    `LOCUS       ${locus} ${construct.length} bp    DNA     linear   SYN`,
    `DEFINITION  CodonForge design ${candidate.candidate_id}.`,
    `ACCESSION   ${candidate.candidate_id}`,
    "VERSION     1",
    "KEYWORDS    codonforge.",
    "SOURCE      synthetic construct",
    "  ORGANISM  synthetic construct",
    "FEATURES             Location/Qualifiers",
    `     source          1..${construct.length}`,
    '                     /mol_type="other DNA"',
    `     CDS             ${five.length + 1}..${five.length + cds.length}`,
    "                     /codon_start=1",
    `                     /transl_table=${sequenceDocument.translation_table}`,
    `                     /note="origin=${candidate.provenance.origin}"`,
    `                     /note="front=${front.front_id}"`,
  ];
  for (let index = 0; index < front.objective_order.length; index += 1) {
    const term = front.objective_order[index];
    const value = formatNumber(candidate.objective_vector[index]);
    const direction = front.objective_sense[index];
    lines.push(`                     /note="${term}=${value} (${direction})"`);
  }
  lines.push("ORIGIN");
  for (let offset = 0; offset < construct.length; offset += 60) {
    const chunk = construct.slice(offset, offset + 60).toLowerCase();
    const blocks = [];
    for (let index = 0; index < chunk.length; index += 10) blocks.push(chunk.slice(index, index + 10));
    lines.push(`${String(offset + 1).padStart(9, " ")} ${blocks.join(" ")}`);
  }
  lines.push("//");
  return `${lines.join("\n")}\n`;
}
