# wsabi
wsabi is a serverless suite for locus-centric exploration of SNVs and CpG methylation. It integrates multi-platform WGS (Illumina, PacBio, ONT) to support genome engineering, cell-line authentication, and monitoring of genetic/epigenetic stability in human iPSC models.


function findNearestSNP(locus) {
    // Filters variants for the correct chromosome and sorts them
    const candidates = STATE.snpListFiltered.filter(s => s.chr === locus.chr).sort((a, b) => a.pos - b.pos);
    if (!candidates.length) return { list: [], index: -1 };

    // Performs Binary Search to find the closest coordinate efficiently
    let lo = 0, hi = candidates.length - 1;
    while (lo <= hi) {
        const mid = (lo + hi) >> 1;
        if (candidates[mid].pos < locus.pos) lo = mid + 1;
        else hi = mid - 1;
    }

    // Resolves whether the 'left' or 'right' hit is geographically closer
    let idx = (hi >= 0 && (Math.abs(candidates[hi].pos - locus.pos) < Math.abs(candidates[lo]?.pos - locus.pos))) ? hi : lo;
    return { list: candidates, index: idx };
}

async function fetchMethContext() {
    const m = METH.data[METH.index];
    // Fetches sequence from API to determine context at the methylated site
    const methylPos = m.pos + 1; 
    const seq = await fetchSequence(m.chr, methylPos - 1, methylPos + 1);

    const methylatedBase = seq[1]; // The middle base at the target position
    let strand = (methylatedBase === 'C') ? '+' : (methylatedBase === 'G') ? '-' : '?';

    // Updates state for visualization in the app
    METH.currentDinuc = { methylatedBase, strand, methylPos };
    METH.target = { chr: m.chr, pos: methylPos };
}
