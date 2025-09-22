function populateKnownVirusesDashboard(data) {
const container = document.getElementById('known-viruses-container');
const noKnownVirusesElement = document.getElementById('no-known-viruses');

// Filter viral hits (≥90% identity, ≥70% coverage in BLASTx)
const knownViruses = data.Viral_Hits.filter(hit => 
hit.BLASTx_Ident >= 90 && hit.BLASTx_Cover >= 70
);

// show or hide the "no known viruses" message
if (knownViruses.length === 0) {
noKnownVirusesElement.style.display = 'block';
return;
} else {
noKnownVirusesElement.style.display = 'none';
}

// card for each known virus
knownViruses.forEach(virus => {
const card = document.createElement('div');
card.className = 'known-virus-card';

// Virus icon 
const virusIcon = `
    <svg viewBox="0 0 24 24" width="24" height="24">
        <circle cx="12" cy="12" r="4" fill="currentColor"/>
        <line x1="12" y1="2" x2="12" y2="6" stroke="currentColor" stroke-width="2"/>
        <line x1="12" y1="18" x2="12" y2="22" stroke="currentColor" stroke-width="2"/>
        <line x1="2" y1="12" x2="6" y2="12" stroke="currentColor" stroke-width="2"/>
        <line x1="18" y1="12" x2="22" y2="12" stroke="currentColor" stroke-width="2"/>
        <line x1="4" y1="4" x2="7" y2="7" stroke="currentColor" stroke-width="2"/>
        <line x1="17" y1="17" x2="20" y2="20" stroke="currentColor" stroke-width="2"/>
        <line x1="4" y1="20" x2="7" y2="17" stroke="currentColor" stroke-width="2"/>
        <line x1="17" y1="7" x2="20" y2="4" stroke="currentColor" stroke-width="2"/>
    </svg>
`;

// Format the organism name for display
const organismName = virus.BLASTx_Organism_Name || virus.ScientificName || "Unknown virus";

// card content
card.innerHTML = `
    <div class="virus-header">
        <div class="virus-icon">${virusIcon}</div>
        <div class="virus-name">${organismName}</div>
    </div>
    
    <div class="virus-metrics">
        <div class="virus-metric">
            <div class="metric-value">${virus.BLASTx_Ident}%</div>
            <div class="metric-label">Identity</div>
        </div>
        <div class="virus-metric">
            <div class="metric-value">${virus.BLASTx_Cover}%</div>
            <div class="metric-label">Coverage</div>
        </div>
        <div class="virus-metric">
            <div class="metric-value">${Math.trunc(virus.BLASTx_Qlength)} nt</div>
            <div class="metric-label">Length</div>
        </div>
    </div>
    
    <div class="virus-taxonomy">
        ${virus.Family ? `<span>Family: ${virus.Family}</span>` : ''}
        ${virus.Genome ? `<span>Genome: ${virus.Genome}</span>` : ''}
    </div>
    
    <a class="contig-link" onclick="scrollToContig('${virus.QueryID}')">
        View details (${virus.QueryID})
    </a>
`;

container.appendChild(card);
});

// Update the stats in the main dashboard
// document.getElementById('viral-sequences').textContent = data.Viral_Hits.length;
// document.getElementById('known-viruses-count').textContent = knownViruses.length;
}

// Helper function to format E-values nicely
function formatEvalue(evalue) {
if (evalue === null || evalue === undefined) return 'N/A';

// If in scientific notation as a string, return it
if (typeof evalue === 'string' && evalue.includes('e')) return evalue;

// Convert to number if it's a string
const evalNum = typeof evalue === 'string' ? parseFloat(evalue) : evalue;

// Format small numbers in scientific notation
if (evalNum < 0.001) {
return evalNum.toExponential(1);
}

return evalNum.toString();
}

// Function to scroll to a specific contig when clicking on a virus card
function scrollToContig(contigId) {
document.querySelector(`.visualization-wrapper[data-contig="${contigId}"]`)
.scrollIntoView({ behavior: 'smooth' });
}


function toggleIndex() {
const index = document.getElementById('contig-index');
const container = document.querySelector('.page-container');
index.classList.toggle('collapsed');
container.classList.toggle('index-collapsed');
}


function createIndex(data) {
const indexList = document.getElementById('index-list');
indexList.innerHTML = ''; // Clear existing items

// Group contigs by sample name
const sampleGroups = {};

// Group viral hits by sample name
data.Viral_Hits.forEach(hit => {
const sampleName = hit.Sample_name || 'Unknown Sample';
if (!sampleGroups[sampleName]) {
    sampleGroups[sampleName] = [];
}
sampleGroups[sampleName].push(hit);
});

// Create expandable sections for each sample
Object.keys(sampleGroups).sort().forEach(sampleName => {
// Create sample group container
const sampleGroup = document.createElement('div');
sampleGroup.className = 'sample-group';

// Create sample header
const sampleHeader = document.createElement('div');
sampleHeader.className = 'sample-header';

// Create toggle icon
const toggleIcon = document.createElement('span');
toggleIcon.className = 'toggle-icon';
toggleIcon.innerHTML = '▼'; // Down arrow for expanded state

// Create sample name element
const sampleNameElement = document.createElement('span');
sampleNameElement.className = 'sample-name';
sampleNameElement.textContent = sampleName;

// Add count badge
const countBadge = document.createElement('span');
countBadge.className = 'count-badge';
countBadge.textContent = sampleGroups[sampleName].length;

// Assemble sample header
sampleHeader.appendChild(toggleIcon);
sampleHeader.appendChild(sampleNameElement);
sampleHeader.appendChild(countBadge);

// Create contig list container for this sample
const contigList = document.createElement('ul');
contigList.className = 'contig-list';

// Add all contigs for this sample
sampleGroups[sampleName].forEach(hit => {
    const li = document.createElement('li');
    li.className = 'contig-item';
    
    // Create a wrapper for the QueryID to allow styling
    const queryIdSpan = document.createElement('span');
    queryIdSpan.className = 'query-id';
    queryIdSpan.textContent = hit.QueryID;
    
    // Add indicators for known viruses (high confidence hits)
    if (hit.vq_score >= 90) {
        const indicator = document.createElement('span');
        indicator.className = 'known-virus-indicator';
        indicator.title = 'Known virus (≥90% identity, ≥70% coverage)';
        indicator.innerHTML = '★'; // Star icon for known viruses
        li.appendChild(indicator);
    }

    if (hit.vq_score <= 20) {
        const indicator = document.createElement('span');
        indicator.className = 'no-virus-indicator';
        // indicator.title = 'Known virus (≥90% identity, ≥70% coverage)';
        indicator.innerHTML = '✘'; // X icon for no viruses
        li.appendChild(indicator);
    }

    
    li.appendChild(queryIdSpan);
    
    // Add click event to scroll to the contig visualization
    li.onclick = () => {
        document.querySelector(`.visualization-wrapper[data-contig="${hit.QueryID}"]`)
            .scrollIntoView({ behavior: 'smooth' });
    };
    
    contigList.appendChild(li);
});

// Toggle sample group expansion when header is clicked
sampleHeader.onclick = (e) => {
    // Don't trigger if clicking on a child element that has its own click handler
    if (e.target !== sampleHeader && e.target !== toggleIcon && e.target !== sampleNameElement) {
        return;
    }
    
    sampleGroup.classList.toggle('collapsed');
    toggleIcon.innerHTML = sampleGroup.classList.contains('collapsed') ? '►' : '▼';
};

// Assemble sample group
sampleGroup.appendChild(sampleHeader);
sampleGroup.appendChild(contigList);
indexList.appendChild(sampleGroup);
});
}

// Add the necessary CSS for the new index structure
const style = document.createElement('style');
style.textContent = `
.sample-group {
margin-bottom: 10px;
}

.sample-header {
display: flex;
align-items: center;
padding: 8px 12px;
background-color: rgba(0, 0, 0, 0.05);
border-radius: 4px;
cursor: pointer;
font-weight: 600;
color: var(--text-primary);
transition: background-color 0.2s;
}

.sample-header:hover {
background-color: rgba(0, 0, 0, 0.1);
}

.toggle-icon {
margin-right: 8px;
font-size: 10px;
transition: transform 0.2s;
}

.sample-name {
flex: 1;
}

.count-badge {
background-color: var(--primary-color);
color: white;
border-radius: 12px;
padding: 2px 8px;
font-size: 12px;
font-weight: 500;
}

.contig-list {
list-style: none;
padding-left: 20px;
margin: 5px 0;
max-height: 60%;
overflow-y: auto;
}

.sample-group.collapsed .contig-list {
display: none;
}

.contig-item {
padding: 6px 10px;
border-radius: 3px;
cursor: pointer;
transition: background-color 0.2s;
display: flex;
align-items: center;
}

.contig-item:hover {
background-color: rgba(0, 0, 0, 0.05);
}

.query-id {
flex: 1;
}

.known-virus-indicator {
color: #f39c12;
margin-right: 6px;
font-size: 14px;
}

.no-virus-indicator {
color: #f32c12;
margin-right: 6px;
font-size: 14px;
}
`;
document.head.appendChild(style);


function createVisualizations(data) {
// container for all visualizations
const container = d3.select("#visualization-container");
container.selectAll("*").remove();

// styling constants
const styles = {
    trackHeight: 130, // accommodate Pfam domains
    orfHeight: 20,
    pfamHeight: 20,
    pfamRadius: 10,
    margin: { top: 50, right: 50, bottom: 100, left: 50 },
    colors: {
        complete: "#3498db",
        "5-prime-partial": "#e74c3c",
        "3-prime-partial": "#f1c40f",
        genomeLine: "#2c3e50",
        axis: "#7f8c8d"
    },
    fonts: {
        primary: "system-ui, -apple-system, sans-serif",
        size: {
            title: "19px",
            subtitle: "15px",
            label: "15px",
            tooltip: "15px"
        }
    },
    pfamColors: [
        "#ff6b6b", "#58b368", "#bc6ff1", "#ffa048",
        "#5191d1", "#ffcc29", "#f06595", "#38b2ac"
    ]
};

    // process each viral hit
    data.Viral_Hits.forEach((viralHit, index) => {
    const contigId = viralHit.QueryID;
    const contigLength = Math.trunc(viralHit.BLASTx_Qlength || viralHit.BLASTn_Qlength);
    const orfs = data.ORF_Data[contigId]
    const organism = viralHit.BLASTx_Organism_Name
    const fullSeq = viralHit.FullSeq
    const vqScore = viralHit.vq_score;


    // individual visualization container
    const visualizationDiv = container.append("div")
        .attr("class", "visualization-wrapper")
        .attr("data-contig", contigId)
        .style("margin-bottom", "40px")
        .style("border-bottom", "1px solid #ecf0f1")
        .style("padding-bottom", "20px");

    // dimensions
    const width = Math.min(window.innerWidth - 40, 1300);
    const height = styles.trackHeight + styles.margin.top + styles.margin.bottom;
    const innerWidth = width - styles.margin.left - styles.margin.right;

    // SVG for this contig
    const svg = visualizationDiv.append("svg")
        .attr("width", width)
        .attr("height", height);

    // main group
    const mainGroup = svg.append("g")
        .attr("transform", `translate(${styles.margin.left}, ${styles.margin.top})`);

    // group for title and button
    const titleGroup = mainGroup.append("g")
        .attr("class", "title-group");

    // title
    const titleText = titleGroup.append("text")
        .attr("class", "visualization-title")
        .attr("x", 0)
        .attr("y", -25)
        .attr("text-anchor", "start")
        .style("font-family", styles.fonts.primary)
        .style("font-size", styles.fonts.size.title)
        .style("font-weight", "600")
        .text(`${contigId}: ${contigLength} nt - ${organism}`);

    // add subtitle
    const subtitleText = titleGroup.append("text")
        .attr("class", "visualization-title")
        .attr("x", 0)
        .attr("y", 0)
        .attr("text-anchor", "start")
        .style("font-family", styles.fonts.primary)
        .style("font-size", styles.fonts.size.subtitle)
        .style("font-weight", "300")
        .text(`Score: ${vqScore}`);

    // bounding box of the title to position the button
    const titleBBox = titleText.node().getBBox();
    const subtitleBBox = subtitleText.node().getBBox();
    const buttonXPosition = titleBBox.width + 10;

    // copy button
    const copyButton = titleGroup.append("g")
        .attr("class", "copy-fasta-button")
        .style("cursor", "pointer")
        .on("click", function() {
            const fastaHeader = `>${contigId}`;
            const fastaSequence = fullSeq;
            const fastaString = `${fastaHeader}\n${fastaSequence}`;

            navigator.clipboard.writeText(fastaString)
                .then(() => {
                    // Feedback visual (opcional)
                    d3.select(this).select("rect")
                        .transition()
                        .duration(200)
                        .style("fill", "#5cb85c")
                        .transition()
                        .duration(1000)
                        .style("fill", getComputedStyle(this).getPropertyValue('--card-background'));
                    d3.select(this).select("text")
                        .transition()
                        .duration(200)
                        .style("fill", "white")
                        .transition()
                        .duration(1000)
                        .style("fill", getComputedStyle(this).getPropertyValue('--primary-color'));
                })
                .catch(err => {
                    console.error('Falha ao copiar: ', err);
                    d3.select(this).select("rect")
                        .style("fill", "#d9534f");
                    d3.select(this).select("text")
                        .style("fill", "white");
                });
        });

    const buttonPaddingHorizontal = 12;
    const buttonPaddingVertical = 6;
    const fontSize = "1rem";
    const borderRadius = 4;
    const textColor = getComputedStyle(document.documentElement).getPropertyValue('--primary-color');
    const backgroundColor = getComputedStyle(document.documentElement).getPropertyValue('--card-background');
    const borderColor = getComputedStyle(document.documentElement).getPropertyValue('--border-color');
    const hoverBackgroundColor = getComputedStyle(document.documentElement).getPropertyValue('--accent-color');
    const hoverTextColor = "white";
    const hoverBorderColor = getComputedStyle(document.documentElement).getPropertyValue('--accent-color');

    const buttonText = "FASTA to clipboard";
    const textElement = copyButton.append("text")
        .attr("x", 0)
        .attr("y", 0)
        .attr("dy", "0.35em")
        .style("text-anchor", "middle")
        .style("font-family", styles.fonts.primary)
        .style("font-size", fontSize)
        .style("fill", textColor)
        .text(buttonText);

    const textWidth = textElement.node().getBBox().width;
    const buttonWidth = textWidth + 2 * buttonPaddingHorizontal;
    const buttonHeight = parseFloat(fontSize) * 1.1 + 5  * buttonPaddingVertical;

    const buttonRect = copyButton.insert("rect", ":first-child")
        .attr("width", buttonWidth)
        .attr("height", buttonHeight)
        .attr("x", -buttonWidth / 2)
        .attr("y", -buttonHeight / 2 + 2)
        .style("fill", backgroundColor)
        .style("stroke", borderColor)
        .style("stroke-width", 1)
        .attr("rx", borderRadius)
        .attr("ry", borderRadius);

    // vertical button adjust
    copyButton.attr("transform", `translate(${innerWidth - buttonWidth - styles.margin.right + 20}, -35)`); // margin adjust

    // adding D3 styles
    copyButton.on("mouseover", function() {
        d3.select(this).select("rect")
            .transition()
            .duration(200)
            .style("fill", hoverBackgroundColor)
            .style("stroke", hoverBorderColor);
        d3.select(this).select("text")
            .transition()
            .duration(200)
            .style("fill", hoverTextColor);
    })
    .on("mouseout", function() {
        d3.select(this).select("rect")
            .transition()
            .duration(200)
            .style("fill", backgroundColor)
            .style("stroke", borderColor);
        d3.select(this).select("text")
            .transition()
            .duration(200)
            .style("fill", textColor);
    });


    // scales
    const xScale = d3.scaleLinear()
        .domain([0, parseInt(contigLength)])
        .range([0, innerWidth]);

    // track group
    const trackGroup = mainGroup.append("g");

    // genome line with gradient
    const gradientId = `genomeLineGradient-${contigId}`;
    const gradient = svg.append("defs")
        .append("linearGradient")
        .attr("id", gradientId)
        .attr("x1", "0%")
        .attr("x2", "100%");
    gradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", "#2c3e50");
    gradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", "#34495e");

    trackGroup.append("line")
        .attr("class", "genome-line")
        .attr("x1", 0)
        .attr("x2", innerWidth)
        .attr("y1", styles.trackHeight / 2)
        .attr("y2", styles.trackHeight / 2)
        .style("stroke", "#2c3e50")
        .style("stroke-width", 2);

    const xAxis = d3.axisBottom(xScale)
        .ticks(10)
        .tickFormat(d => `${d.toLocaleString()}bp`);

    trackGroup.append("g")
        .attr("class", "scale-axis")
        .attr("transform", `translate(0, ${styles.trackHeight + 10})`)
        .call(xAxis)
        .style("font-family", styles.fonts.primary)
        .style("font-size", styles.fonts.size.label);



    // BLAST and Taxonomy information panel 
    const infoPanel = visualizationDiv.append("div")
        .attr("class", "info-panel")
        .style("margin-bottom", "20px")
        .style("padding", "15px")
        .style("border", "1px solid #ecf0f1")
        .style("border-radius", "4px")
        .style("background-color", "#f9f9f9");

    // two-column layout for BLAST info
    const blastInfo = infoPanel.append("div")
        .style("display", "grid")
        .style("grid-template-columns", "1fr 1fr")
        .style("gap", "15px")
        .style("margin-bottom", "15px");

    // BLASTx column
    blastInfo.append("div")
        .html(`
            <strong style="color: #2c3e50; font-size: 14px">BLASTx Results</strong>
            <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
            ${viralHit.BLASTx_Subject_Title ? `
                <strong>Subject:</strong> ${viralHit.BLASTx_Subject_Title}<br>
                <strong>Organism:</strong> ${viralHit.BLASTx_Organism_Name}<br>
                <strong>Coverage:</strong> ${viralHit.BLASTx_Cover}%<br>
                <strong>Identity:</strong> ${viralHit.BLASTx_Ident}%<br>
                <strong>E-value:</strong> ${viralHit.BLASTx_evalue}<br>
                <strong>Subj. Length:</strong> ${Math.trunc(viralHit.BLASTx_Slength)} aa
            ` : 'No BLASTx hits'}
        `);

    // BLASTn column
    blastInfo.append("div")
        .html(`
            <strong style="color: #2c3e50; font-size: 14px">BLASTn Results</strong>
            <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
            ${viralHit.BLASTn_Subject_Title ? `
                <strong>Subject:</strong> ${viralHit.BLASTn_Subject_Title}<br>
                <strong>Coverage:</strong> ${viralHit.BLASTn_Cover}%<br>
                <strong>Identity:</strong> ${viralHit.BLASTn_Ident}%<br>
                <strong>E-value:</strong> ${viralHit.BLASTn_evalue}<br>
                <strong>Subj. Length:</strong> ${Math.trunc(viralHit.BLASTn_Slength)} nt
            ` : 'No BLASTn hits'}
        `);

    // taxonomy section
    infoPanel.append("div")
        .style("grid-column", "span 2")
        .style("margin-bottom", "15px")
        .style("padding", "15px")
        .style("border", "1px solid #ecf0f1")
        .style("border-radius", "5px")
        .style("background-color", "#f9f9f9")
        .html(`
            <strong style="color: #2c3e50; font-size: 16px; display: block; margin-bottom: 10px;">Taxonomy - BLASTx Hit</strong>
            <hr style="border: 1px solid #ddd; margin: 10px 0;">
            <div style="margin-bottom: 8px;">
                <strong style="color: #555; font-weight: 600; margin-right: 5px;">Scientific Name:</strong> ${viralHit.ScientificName}
            </div>
            <div style="margin-bottom: 8px;">
                <strong style="color: #555; font-weight: 600; margin-right: 5px;">Species:</strong> <span style="color: ${viralHit.Species ? '#333' : '#999'}; font-style: ${viralHit.Species ? 'normal' : 'italic'};">${viralHit.Species || 'N/A'}</span>
            </div>
            <div style="margin-bottom: 8px;">
                <strong style="color: #555; font-weight: 600; margin-right: 5px;">Rank:</strong> <span style="color: ${viralHit.NoRank ? '#333' : '#999'}; font-style: ${viralHit.NoRank ? 'normal' : 'italic'};">${viralHit.NoRank || 'N/A'}</span>
            </div>
            <br>
            <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 10px;">                           
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">Kingdom</strong>
                    <span style="color: ${viralHit.Kingdom ? '#333' : '#999'}; font-style: ${viralHit.Kingdom ? 'normal' : 'italic'}; display: block;">${viralHit.Kingdom || 'N/A'}</span>
                </div>
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">Phylum</strong>
                    <span style="color: ${viralHit.Phylum ? '#333' : '#999'}; font-style: ${viralHit.Phylum ? 'normal' : 'italic'}; display: block;">${viralHit.Phylum || 'N/A'}</span>
                </div>
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">Class</strong>
                    <span style="color: ${viralHit.Class ? '#333' : '#999'}; font-style: ${viralHit.Class ? 'normal' : 'italic'}; display: block;">${viralHit.Class || 'N/A'}</span>
                </div>
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">Order</strong>
                    <span style="color: ${viralHit.Order ? '#333' : '#999'}; font-style: ${viralHit.Order ? 'normal' : 'italic'}; display: block;">${viralHit.Order || 'N/A'}</span>
                </div>
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">Family</strong>
                    <span style="color: ${viralHit.Family ? '#333' : '#999'}; font-style: ${viralHit.Family ? 'normal' : 'italic'}; display: block;">${viralHit.Family || 'N/A'}</span>
                </div>
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">Genome</strong>
                    <span style="color: ${viralHit.Genome ? '#333' : '#999'}; font-style: ${viralHit.Genome ? 'normal' : 'italic'}; display: block;">${viralHit.Genome || 'N/A'}</span>
                </div>
                <div>
                    <strong style="color: #555; font-weight: 600; display: block;">TaxID</strong>
                    <span style="color: ${viralHit.TaxId ? '#333' : '#999'}; font-style: ${viralHit.TaxId ? 'normal' : 'italic'}; display: block;">${viralHit.TaxId ? Math.trunc(viralHit.TaxId) : 'N/A'}</span>
                </div>                         
            </div>
        `);

    // Pfam domain summary panel
    const hmmHits = data.HMM_hits.filter(hit => hit.QueryID === contigId);
    if (hmmHits.length > 0) {
        const pfamPanel = infoPanel.append("div")
            .style("margin-top", "15px")
            .style("border-top", "1px solid #ecf0f1")
            .style("padding-top", "15px");

        pfamPanel.append("div")
            .html(`
                <strong style="color: #2c3e50; font-size: 14px">Pfam Conserved Regions</strong>
                <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
            `);

        const pfamList = pfamPanel.append("div")
            .style("display", "grid")
            .style("grid-template-columns", "repeat(auto-fill, minmax(300px, 1fr))")
            .style("gap", "10px");

        // unique Pfam domains across all ORFs in contig
        const pfamDomains = [];
        hmmHits.forEach(hit => {
            // extract keys that match Pfam_Accession
            const pfamKeys = Object.keys(hit).filter(key => key.match(/^Pfam_Accession(_\d+)?$/));

            // Pfam domain information
            pfamKeys.forEach(accKey => {
                const suffix = accKey.replace('Pfam_Accession', '');
                if (hit[accKey] && hit[`Pfam_Description${suffix}`]) {
                    pfamDomains.push({
                        accession: hit[accKey],
                        description: hit[`Pfam_Description${suffix}`],
                        type: hit[`Pfam_Type${suffix}`] || 'Unknown',
                        orf: hit.Query_name
                    });
                }
            });
        });

        // remove duplicates and show unique Pfam domain
        const uniquePfamDomains = Array.from(new Map(
            pfamDomains.map(domain => [domain.accession, domain])
        ).values());

        uniquePfamDomains.forEach((domain, i) => {
            const colorIndex = i % styles.pfamColors.length;
            const domainColor = styles.pfamColors[colorIndex];

            pfamList.append("div")
                .style("padding", "5px")
                .style("border-left", `4px solid ${domainColor}`)
                .style("background-color", "#fff")
                .html(`
                    <strong>${domain.accession}</strong>: ${domain.description}<br>
                    <small>Type: ${domain.type}</small>
                `);
        });
    }

    if (viralHit.AI_summary) {
        const aiSummaryPanel = infoPanel.append("div")
            .style("margin-top", "15px")
            .style("border-top", "1px solid #ecf0f1")
            .style("padding-top", "15px");

        aiSummaryPanel.append("div")
            .html(`
                <strong style="color: #2c3e50; font-size: 16px">AI Summary</strong>
                <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
            `);

        const aiSummaryContent = aiSummaryPanel.append("div")
            .style("padding", "10px")
            .style("background-color", "#fff")
            .style("border-left", "4px solid #9b59b6")
            .style("font-family", styles.fonts.primary)
            .style("line-height", "1.5");

        // Convert Markdown format to HTML
        const markdownText = viralHit.AI_summary;
        
        // Simple markdown conversion for common elements
        const htmlText = markdownText
            .replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>') // Bold text
            .replace(/\*(.*?)\*/g, '<em>$1</em>') // Italic text
            .replace(/^\* (.*?)$/gm, '• $1<br>') // Bullet points
            .replace(/\n\n/g, '<br><br>'); // Line breaks
        
        aiSummaryContent.html(htmlText);
    }

    // tooltip (one per visualization)
    const tooltip = d3.select("body").append("div")
        .attr("class", `tooltip-${contigId}`)
        .style("position", "absolute")
        .style("display", "none")
        .style("background", "white")
        .style("border", "1px solid #ddd")
        .style("border-radius", "4px")
        .style("padding", "10px")
        .style("box-shadow", "2px 2px 6px rgba(0,0,0,0.2)")
        .style("font-family", styles.fonts.primary)
        .style("font-size", styles.fonts.size.tooltip)
        .style("pointer-events", "auto");

    // track tooltip state
    let tooltipTimeout = null;
    let activeTooltip = null;
    let isTooltipPinned = false;

    // create tooltip content for ORFs
    const createTooltipContent = (orf, hmmHit) => {
        const sequence = orf.sequence || '';
        const truncatedSeq = sequence.substring(0, 8) + "..";
        const buttonId = `copy-button-${orf.Query_name.replace(/[^a-zA-Z0-9]/g, '-')}`;
        const feedbackId = `feedback-${orf.Query_name.replace(/[^a-zA-Z0-9]/g, '-')}`;


        // position string based on strand
        const positionString = orf.strand === '+' ? 
            `${orf.start}-${orf.end}` : 
            `${orf.end}-${orf.start}`;

        return `
            <div style="font-family: ${styles.fonts.primary}">
                <strong style="color: #2c3e50">${orf.Query_name}</strong><br>
                <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                <strong>Position (nt):</strong> ${positionString}<br>
                <strong>Type:</strong> ${orf.type}<br>
                <strong>Length (nt):</strong> ${orf.length}<br>
                <strong>Length (aa):</strong> ${orf.length / 3}<br>
                <strong>Frame:</strong> ${orf.frame}<br>
                <strong>Strand:</strong> ${orf.strand}<br>
                <strong>Codons:</strong> ${orf.start_codon} → ${orf.stop_codon}<br>
                <strong>Sequence (aa):</strong>
                <span style="font-family: monospace">${truncatedSeq}</span>
                <button
                    id="${buttonId}"
                    class="copy-button"
                    data-sequence="${sequence}"
                    onclick="handleCopyClick(this)">
                    Copy
                </button>
                <span id="${feedbackId}" class="copy-feedback">
                    Copied!
                </span>
                ${hmmHit ? `
                    <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                    <strong style="color: #2c3e50">HMM Hits:</strong><br>
                    ${['RVDB', 'Vfam', 'EggNOG'].map(db =>
                        hmmHit[`${db}_TargetID`] ?
                            `<strong>${db}:</strong> ${hmmHit[`${db}_TargetID`]}
                            ${hmmHit[`${db}_Description`] ? ` - ${hmmHit[`${db}_Description`]}` : ''}<br>`
                        : ''
                    ).join('')}
                ` : ''}
            </div>
        `;
    };

    // tooltip content for Pfam domains
    const createPfamTooltipContent = (pfamData) => {
        return `
            <div class="pfam-tooltip" style="font-family: ${styles.fonts.primary}">
                <strong style="color: #2c3e50">${pfamData.Pfam_Accession}</strong><br>
                <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                <strong>Description:</strong> ${pfamData.Pfam_Description}<br>
                <strong>Type:</strong> ${pfamData.Pfam_Type}<br>
                <strong>Position (aa):</strong> ${pfamData.Pfam_Start}-${pfamData.Pfam_End}<br>
                <strong>Length (aa):</strong> ${pfamData.Pfam_length}<br>
                <strong>Score:</strong> ${parseFloat(pfamData.Pfam_Score).toFixed(2)}<br>
                ${pfamData.Pfam_Info ? `
                    <div class="pfam-info">
                        ${pfamData.Pfam_Info.replace(/<\/?p>/g, '')}
                    </div>
                ` : ''}
            </div>
        `;
    };

    // handle copy button click
    window.handleCopyClick = function(button) {
        const sequence = button.getAttribute('data-sequence');
        const feedbackId = button.id.replace('copy-button-', 'feedback-');
        const feedbackElement = document.getElementById(feedbackId);

        navigator.clipboard.writeText(sequence)
            .then(() => {
                // feedback
                button.style.display = 'none';
                feedbackElement.style.display = 'inline';
                feedbackElement.classList.add('show');

                // reset after 2 seconds
                setTimeout(() => {
                    button.style.display = 'inline';
                    feedbackElement.classList.remove('show');
                    setTimeout(() => {
                        feedbackElement.style.display = 'none';
                    }, 300);
                }, 2000);
            })
            .catch(err => {
                console.error('Failed to copy:', err);
                feedbackElement.textContent = 'Error copying!';
                feedbackElement.style.color = '#e74c3c';
                feedbackElement.style.display = 'inline';
                feedbackElement.classList.add('show');
            });
    };

    // show tooltip
    const showTooltip = (event, content) => {
        tooltip
            .style("display", "block")
            .style("left", (event.pageX + 15) + "px")
            .style("top", (event.pageY - 10) + "px")
            .html(content);
    };

    // hide tooltip
    const hideTooltip = () => {
        if (!isTooltipPinned) {
            tooltip.style("display", "none");
            activeTooltip = null;
        }
    };

    // process ORFs
    orfs.forEach(orf => {
        const orfGroup = trackGroup.append("g");
        const yOffset = orf.frame > 0 ? -0 : 0;
        const orfWidth = xScale(orf.end) - xScale(orf.start);
        const arrowSize = Math.min(15, orfWidth / 5);

        // arrow-shaped ORF path
        const orfPath = orf.strand === "+" ?
            `M ${xScale(orf.start)} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset}
            L ${xScale(orf.end) - arrowSize} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset}
            L ${xScale(orf.end)} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight / 2}
            L ${xScale(orf.end) - arrowSize} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}
            L ${xScale(orf.start)} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}
            Z` :
            `M ${xScale(orf.start)} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight / 2}
            L ${xScale(orf.start) + arrowSize} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset}
            L ${xScale(orf.end)} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset}
            L ${xScale(orf.end)} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}
            L ${xScale(orf.start) + arrowSize} ${(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}
            Z`;

        // ORF shape with interaction
        orfGroup.append("path")
            .attr("d", orfPath)
            .style("fill", styles.colors[orf.type])
            .style("filter", "none")
            .style("stroke", "#333")
            .style("cursor", "pointer")
            .on("mouseover", (event) => {
                const hmmHit = data.HMM_hits.find(hit => hit.Query_name === orf.Query_name);
                isTooltipPinned = false;
                showTooltip(event, createTooltipContent(orf, hmmHit));
            })
            .on("mousemove", (event) => {
                if (!isTooltipPinned) {
                    tooltip
                        .style("left", (event.pageX + 15) + "px")
                        .style("top", (event.pageY - 10) + "px");
                }
            })
            .on("mouseout", () => {
                tooltipTimeout = setTimeout(hideTooltip, 3000);
            })
            .on("click", (event) => {
                const hmmHit = data.HMM_hits.find(hit => hit.Query_name === orf.Query_name);
                isTooltipPinned = !isTooltipPinned;
                if (isTooltipPinned) {
                    showTooltip(event, createTooltipContent(orf, hmmHit));
                    if (tooltipTimeout) {
                        clearTimeout(tooltipTimeout);
                    }
                } else {
                    hideTooltip();
                }
                event.stopPropagation();
            });

        // corresponding HMM hit and extract Pfam domains
        const hmmHit = data.HMM_hits.find(hit => hit.Query_name === orf.Query_name);
        if (hmmHit) {
            // all Pfam domains from the HMM hit
            const pfamDomains = [];

            // keys that match Pfam_Accession pattern
            const pfamKeys = Object.keys(hmmHit).filter(key => key.match(/^Pfam_Accession(_\d+)?$/));

            // domain object with all related information
            pfamKeys.forEach(accKey => {
                // suffix (empty string or _2, _3, etc.)
                const suffix = accKey.replace('Pfam_Accession', '');

                // process if there's a start and end position
                if (hmmHit[`Pfam_Start${suffix}`] && hmmHit[`Pfam_End${suffix}`]) {
                    pfamDomains.push({
                        Pfam_Accession: hmmHit[accKey],
                        Pfam_Description: hmmHit[`Pfam_Description${suffix}`] || '',
                        Pfam_Info: hmmHit[`Pfam_Info${suffix}`] || '',
                        Pfam_Type: hmmHit[`Pfam_Type${suffix}`] || '',
                        Pfam_Score: hmmHit[`Pfam_Score${suffix}`] || '0',
                        Pfam_Start: parseFloat(hmmHit[`Pfam_Start${suffix}`]),
                        Pfam_End: parseFloat(hmmHit[`Pfam_End${suffix}`]),
                        Pfam_length: parseFloat(hmmHit[`Pfam_length${suffix}`] || '0')
                    });
                }
            });

            // Pfam domains below the ORF
            const pfamGroup = trackGroup.append("g");

            pfamDomains.forEach((domain, i) => {
                const colorIndex = i % styles.pfamColors.length;
                const domainColor = styles.pfamColors[colorIndex];

                // Pfam domain positions from amino acid to nucleotide coordinates
                let domainStartNuc, domainEndNuc;

                // frame offset (0, 1, or 2 depending on the reading frame)
                const frameOffset = Math.abs(orf.frame) - 1; // Convert frame (1-6) to offset (0-2)

                if (orf.strand === "+") {
                    // for positive strand, translate directly from AA to nucl
                    // each amino acid is 3 nucleotides
                    domainStartNuc = orf.start + frameOffset + (domain.Pfam_Start - 1) * 3;
                    domainEndNuc = orf.start + frameOffset + (domain.Pfam_End * 3) - 1; // -1 because Pfam_End is inclusive
                } else {
                    // for strand == -, we count from the end backwards
                    // We need to reverse the amino acid coordinates
                    domainStartNuc = orf.end - frameOffset - (domain.Pfam_End * 3) + 1; // +1 because we're counting backwards
                    domainEndNuc = orf.end - frameOffset - (domain.Pfam_Start - 1) * 3;
                }

                // ensure domain stays within the ORF boundaries
                domainStartNuc = Math.max(orf.start, domainStartNuc);
                domainEndNuc = Math.min(orf.end, domainEndNuc);

                const yPos = (styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight + 15;

                pfamGroup.append("rect")
                    .attr("class", "pfam-domain")
                    .attr("x", xScale(domainStartNuc))
                    .attr("y", yPos)
                    .attr("width", xScale(domainEndNuc) - xScale(domainStartNuc))
                    .attr("height", styles.pfamHeight)
                    .attr("rx", styles.pfamRadius)
                    .attr("ry", styles.pfamRadius)
                    .style("fill", domainColor)
                    .style("filter", "none")
                    .style("cursor", "pointer")
                    .on("mouseover", (event) => {
                        showTooltip(event, createPfamTooltipContent(domain));
                    })
                    .on("mousemove", (event) => {
                        if (!isTooltipPinned) {
                            tooltip
                                .style("left", (event.pageX + 15) + "px")
                                .style("top", (event.pageY - 10) + "px");
                        }
                    })
                    .on("mouseout", () => {
                        tooltipTimeout = setTimeout(hideTooltip, 3000);
                    })
                    .on("click", (event) => {
                        isTooltipPinned = !isTooltipPinned;
                        if (isTooltipPinned) {
                            showTooltip(event, createPfamTooltipContent(domain));
                            if (tooltipTimeout) {
                                clearTimeout(tooltipTimeout);
                            }
                        } else {
                            hideTooltip();
                        }
                        event.stopPropagation();
                    });

                // domain label if space allows
                if (xScale(domainEndNuc) - xScale(domainStartNuc) > 80) {
                    pfamGroup.append("text")
                        .attr("x", xScale(domainStartNuc) + (xScale(domainEndNuc) - xScale(domainStartNuc)) / 2)
                        .attr("y", yPos + styles.pfamHeight / 2)
                        .attr("text-anchor", "middle")
                        .attr("dominant-baseline", "central")
                        .style("font-family", styles.fonts.primary)
                        .style("font-size", "10px")
                        .style("fill", "#ffffff")
                        .style("pointer-events", "none")
                        .text(domain.Pfam_Accession);
                }
            });
        }
    });
});

// hide loading indicator
document.getElementById("loading").style.display = "none";

// add event listener to close tooltip when clicking elsewhere
document.addEventListener("click", (event) => {
    if (!event.target.closest(".tooltip")) {
        d3.selectAll("[class^='tooltip-']").style("display", "none");
        isTooltipPinned = false;
    }
});
}

// load data and initialize visualization
document.addEventListener("DOMContentLoaded", () => {


// data from previous code
const logData = {
number_originalSeqs: 216,
number_filteredSeqs: 216,
number_viralSeqs: 59,
args: {
    input: "seq.fasta",
    outdir: "teste_newPfam",
    numORFs: 665,
    cap3: true,
    cpu: 4,
    blastn: "viral_genomic.fna",
    diamond_blastx: "viralDB.dmnd",
    pfam_hmm: "Pfam-A.hmm"
}
};

// populate the dashboard
function populateDashboard(data) {
// statistic cards
document.getElementById('total-sequences').textContent = data.number_originalSeqs.toLocaleString();
document.getElementById('sequences-500nt').textContent = data.number_filteredSeqs.toLocaleString();
document.getElementById('viral-sequences').textContent = data.number_viralSeqs.toLocaleString();
document.getElementById('num-orfs').textContent = data.args.numORFs.toLocaleString();

// configuration information
document.getElementById('input-file').textContent = data.args.input;
document.getElementById('output-dir').textContent = data.args.outdir;
document.getElementById('cap3').textContent = data.args.cap3 ? "Enabled" : "Disabled";
document.getElementById('cpu-cores').textContent = data.args.cpu;

// database information
document.getElementById('blastn-db').textContent = `BLASTn: ${data.args.blastn}`;
document.getElementById('blastx-db').textContent = `BLASTx: ${data.args.diamond_blastx}`;
document.getElementById('pfam_hmm').textContent = `Pfam: ${data.args.pfam_hmm}`;

// pie chart for sequence distribution
createSequenceDistributionChart(data);
}

// create sequence distribution chart
function createSequenceDistributionChart(data) {
const width = document.getElementById('sequences-chart').clientWidth;
const height = 200;
const radius = Math.min(width, height) / 2 * 0.8;

const svg = d3.select("#sequences-chart")
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("transform", `translate(${width/2}, ${height/2})`);

// data for pie chart
const pieData = [
    { label: "Viral", value: data.number_viralSeqs, color: "#e74c3c" },
    { label: "Non-viral ≥ 500nt", value: data.number_filteredSeqs - data.number_viralSeqs, color: "#3498db" },
    { label: "Contigs < 500nt", value: data.number_originalSeqs - data.number_filteredSeqs, color: "#95a5a6" }
];

// pie chart
const pie = d3.pie()
    .value(d => d.value)
    .sort(null);

const arc = d3.arc()
    .innerRadius(radius * 0.5) // Donut chart
    .outerRadius(radius);

const arcLabels = d3.arc()
    .innerRadius(radius * 0.7)
    .outerRadius(radius * 0.7);

// chart segments
const segments = svg.selectAll("path")
    .data(pie(pieData))
    .enter()
    .append("path")
    .attr("d", arc)
    .attr("fill", d => d.data.color)
    .attr("stroke", "white")
    .style("stroke-width", "2px")
    .style("opacity", 0.8)
    .on("mouseover", function() {
        d3.select(this)
            .style("opacity", 1);
    })
    .on("mouseout", function() {
        d3.select(this)
            .style("opacity", 0.8);
    });

    // background rectangles for percentage labels
    svg.selectAll("rect.label-bg") // Seleciona elementos rect com a classe label-bg (se existirem)
        .data(pie(pieData))
        .enter()
        .append("rect")
        .attr("class", "label-bg") // add class to select style
        .attr("transform", d => `translate(${arcLabels.centroid(d)})`) // arc centroid
        .attr("width", d => { // width based on text
            const percent = (d.data.value / data.number_originalSeqs * 100).toFixed(1);
            const textLength = percent > 0.00000001 ? `${percent}%`.length * 8 : 0; // width estimative 
            return textLength + 10; // add a few padding
        })
        .attr("height", 16) // fix height 
        .attr("x", d => { // adjust x to centralize rectangle
            const percent = (d.data.value / data.number_originalSeqs * 100).toFixed(1);
            const textLength = percent > 0.00000001 ? `${percent}%`.length * 8 : 0;
            return - (textLength + 10) / 2;
        })
        .attr("y", -8) // adjust y to centralize rectangle
        .attr("fill", "white")
        .attr("rx", 5) // radius for round corners
        .attr("ry", 5);

// percentage labels
svg.selectAll("text")
    .data(pie(pieData))
    .enter()
    .append("text")
    .attr("transform", d => `translate(${arcLabels.centroid(d)})`)
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "middle")
    .text(d => {
        const percent = (d.data.value / data.number_originalSeqs * 100).toFixed(1);
        return percent > 0.00000001 ? `${percent}%` : '';
    })
    .style("font-size", "12px")
    .style("fill", "black")
    .style("font-weight", "bold");

// legend
const legend = svg.selectAll(".legend")
    .data(pieData)
    .enter()
    .append("g")
    .attr("class", "legend")
    .attr("transform", (d, i) => `translate(${width/2 - 180}, ${i * 20 - 30})`);

legend.append("rect")
    .attr("width", 12)
    .attr("height", 12)
    .attr("fill", d => d.color);

legend.append("text")
    .attr("x", 20)
    .attr("y", 6)
    .attr("dy", ".35em")
    .style("font-size", "12px")
    .text(d => `${d.label} (${d.value})`);
}





// graph plot
function createTaxonomyVisualization(data) {
// Select the container for the taxonomy visualization
const container = d3.select("#taxonomy-container");
container.selectAll("*").remove();

// Hide the "no known viruses" message if it exists
d3.select("#no-known-viruses").style("display", "none");

// Create a new div for our visualization
const taxonomyDiv = container.append("div")
    .attr("class", "taxonomy-visualization-wrapper")
    .style("width", "97%")
    .style("height", "500px")
    .style("background-color", "#ffffff")
    .style("border-radius", "8px")
    .style("box-shadow", "0 2px 5px rgba(0, 0, 0, 0.05)")
    .style("padding", "20px")
    .style("margin-bottom", "20px");


// SVG container for the visualization
const width = taxonomyDiv.node().getBoundingClientRect().width;
const height = 500;

const svg = taxonomyDiv.append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("viewBox", `0 0 ${width} ${height-210}`)
    .attr("preserveAspectRatio", "xMidYMid meet");

// Main group
const mainGroup = svg.append("g")
    .attr("transform", `translate(${width / 2}, ${height / 3.4})`);

// Enhanced styling constants for more compact visualization
const styles = {
    margin: { top: 100, right: 10, bottom: 10, left: 10 },
    colors: {
        phylum: "#3498db",     // Blue for phylum
        family: "#2ecc71",     // Green for family
        virus: "#e74c3c",      // Red for viruses
        link: "#c0c4c4",       // Light gray for links
        highlight: "#ff7f0e"   // Orange for highlighting
    },
    fonts: {
        primary: "system-ui, -apple-system, sans-serif",
        size: {
            title: "15px",     // Slightly smaller
            label: "11px",     // Slightly smaller
            tooltip: "13px"    // Slightly smaller
        }
    },
    // More subtle color scale for different phyla
    phylumColorScale: d3.scaleOrdinal()
        .range([
            //"#3498db", // Blue
            "#9b59b6", // Purple
            "#e74c3c", // Red
            "#f1c40f", // Yellow
            "#2ecc71", // Green
            "#1abc9c", // Teal
            "#34495e"  // Dark blue
        ]),
    // More subtle color scale for different families
    familyColorScale: d3.scaleOrdinal()
        .range([
            "#2ecc71", // Green
            "#e67e22", // Orange
            "#9b59b6", // Purple
            "#f39c12", // Yellow
            "#16a085", // Teal
            "#d35400", // Dark orange
            "#27ae60", // Dark green
            "#8e44ad"  // Dark purple
        ]),
    // Configuration for force simulation
    force: {
        center: 0.1,         // Strength of center force
        charge: 90,         // Base charge strength
        link: 30,            // Base link distance
        boundary: 0.15,      // Strength of boundary force
        radius: 30         // Boundary radius
    }
};

// Process the data to extract unique phyla, families, and viruses
// This is a placeholder; you'll need to adapt this to your actual data structure
const processedData = processViralData(data);

// Create force simulation with more constraints for compact layout
const simulation = d3.forceSimulation(processedData.nodes)
    .force("link", d3.forceLink(processedData.links).id(d => d.id).distance(d => {
        // Shorter distances for a more compact layout
        if (d.source.type === "phylum" && d.target.type === "family") return 100;
        return 20;
    }))
    .force("charge", d3.forceManyBody().strength(d => {
        // Less repulsion for more compact layout
        if (d.type === "phylum") return -1000;
        if (d.type === "family") return -250;
        return -70;
    }))
    .force("center", d3.forceCenter(0, 0))
    // Add force to keep nodes within bounds
    .force("x", d3.forceX().strength(0.1))
    .force("y", d3.forceY().strength(0.1))
    .force("collision", d3.forceCollide().radius(d => {
        // Smaller collision radiuses for compact layout
        if (d.type === "phylum") return 100;
        if (d.type === "family") return 50;
        return 20;
    }));

// Create links
const link = mainGroup.append("g")
    .attr("class", "links")
    .selectAll("line")
    .data(processedData.links)
    .enter().append("line")
    .attr("stroke", styles.colors.link)
    .attr("stroke-opacity", 0.6)
    .attr("stroke-width", d => {
        // Different stroke widths based on link type
        if (d.source.type === "phylum" && d.target.type === "family") return 2;
        return 1;
    });

// Create nodes
const node = mainGroup.append("g")
    .attr("class", "nodes")
    .selectAll("g")
    .data(processedData.nodes)
    .enter().append("g")
    .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));

// Add circles to nodes - smaller for compactness
node.append("circle")
    .attr("r", d => {
        // Smaller sizes for compact layout
        if (d.type === "phylum") return 20;
        if (d.type === "family") return 12;
        return 5;
    })
    .attr("fill", d => {
        // Different colors based on node type
        if (d.type === "phylum") return styles.phylumColorScale(d.id);
        if (d.type === "family") return styles.familyColorScale(d.id);
        return styles.colors.virus;
    })
    .attr("stroke", "#fff")
    .attr("stroke-width", 1.5);

// Add labels only to phylum and family nodes (omit viral names)
node.filter(d => d.type !== "virus") // Only add text to phylum and family nodes
    .append("text")
    .text(d => d.name)
    .attr("x", 0)
    .attr("y", d => {
        // Position text differently based on node type
        if (d.type === "phylum") return -25;
        if (d.type === "family") return -15;
        return 0;
    })
    .attr("text-anchor", "middle")
    .attr("font-family", styles.fonts.primary)
    .attr("font-size", d => {
        // Different font sizes based on node type
        if (d.type === "phylum") return styles.fonts.size.title;
        return styles.fonts.size.label;
    })
    .attr("fill", "#2c3e50")
    .style("pointer-events", "none");



// Setup tooltip
const tooltip = d3.select("body").append("div")
    .attr("class", "taxonomy-tooltip")
    .style("position", "absolute")
    .style("visibility", "hidden")
    .style("background-color", "#fff")
    .style("border", "1px solid #ddd")
    .style("border-radius", "4px")
    .style("padding", "10px")
    .style("box-shadow", "0 2px 5px rgba(0, 0, 0, 0.1)")
    .style("font-family", styles.fonts.primary)
    .style("font-size", styles.fonts.size.tooltip)
    .style("pointer-events", "none")
    .style("z-index", 1000);

// Add mouseover events for nodes with updated sizing
node.on("mouseover", function(event, d) {
    d3.select(this).select("circle")
        .transition()
        .duration(200)
        .attr("r", d => {
            if (d.type === "phylum") return 25;
            if (d.type === "family") return 16;
            return 8; // Slightly larger for viruses but still compact
        });

    // Show tooltip with information - show viral names here instead of on the graph
    tooltip.style("visibility", "visible")
        .html(() => {
            let content = `<strong>${d.name}</strong><br>Type: ${d.type.charAt(0).toUpperCase() + d.type.slice(1)}`;
            if (d.type === "virus") {
                content += `<br>Family: ${d.family || 'N/A'}`;
                content += `<br>Phylum: ${d.phylum || 'N/A'}`;
                // Add more virus details if available
                if (d.Genome) content += `<br>Genome Type: ${d.Genome}`;
                if (d.taxId) content += `<br>TaxID: ${d.taxId}`;
            } else if (d.type === "family") {
                content += `<br>Phylum: ${d.phylum || 'N/A'}`;
                // Add count of viruses in this family
                const virusCount = processedData.links.filter(link => 
                    link.source.id === d.id || 
                    (typeof link.source === 'string' && link.source === d.id)
                ).length;
                content += `<br>Viruses: ${virusCount}`;
            } else if (d.type === "phylum") {
                // Add count of families in this phylum
                const familyCount = processedData.links.filter(link => 
                    (link.source.type === "phylum" && link.source.id === d.id) || 
                    (typeof link.source === 'string' && link.source === d.id && link.target.type === "family")
                ).length;
                content += `<br>Families: ${familyCount}`;
            }
            return content;
        });

        // Highlight all connected links and nodes
        link.style("opacity", l => 
            (l.source.id === d.id || l.target.id === d.id) ? 1 : 0.1
        )
        .style("stroke-width", l => 
            (l.source.id === d.id || l.target.id === d.id) ? 3 : 1
        )
        .style("stroke", l => 
            (l.source.id === d.id || l.target.id === d.id) ? "#3498db" : styles.colors.link
        );
        
        // Highlight connected nodes
        node.style("opacity", n => 
            (n.id === d.id || 
            processedData.links.some(l => 
                (l.source.id === d.id && l.target.id === n.id) || 
                (l.target.id === d.id && l.source.id === n.id)
            )) ? 1 : 0.3
        );
})
.on("mousemove", function(event) {
    tooltip.style("top", (event.pageY - 10) + "px")
        .style("left", (event.pageX + 10) + "px");
})
.on("mouseout", function(event, d) {
    d3.select(this).select("circle")
        .transition()
        .duration(200)
        .attr("r", d => {
            if (d.type === "phylum") return 20;
            if (d.type === "family") return 12;
            return 5;
        });
    
    tooltip.style("visibility", "hidden");

    // Reset all links and nodes
    link.style("opacity", 0.6)
        .style("stroke-width", d => d.source.type === "phylum" ? 2 : 1)
        .style("stroke", styles.colors.link);
    
    node.style("opacity", 1);
    
    tooltip.style("visibility", "hidden");

});



// constrain the graph in  a rectangular area
const bounds = {
    x: { min: -width/2 + 40, max: width/2 - 40 },
    y: { min: -height/2 + 40, max: height/2 - 20 }
};

simulation.on("tick", () => {
    // Constrain nodes to a rectangular area to prevent flying away
    processedData.nodes.forEach(d => {
        // Apply rectangular constraints
        d.x = Math.max(bounds.x.min, Math.min(bounds.x.max, d.x));
        d.y = Math.max(bounds.y.min, Math.min(bounds.y.max, d.y));
    });
    
    link
        .attr("x1", d => d.source.x)
        .attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x)
        .attr("y2", d => d.target.y);

    node.attr("transform", d => `translate(${d.x},${d.y})`);
});

// Drag functions
function dragstarted(event, d) {
    if (!event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
}

function dragged(event, d) {
    d.fx = event.x;
    d.fy = event.y;
}

function dragended(event, d) {
    if (!event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
}

// Function to process viral data
function processViralData(data) {
    const nodes = [];
    const links = [];
    const phylaMap = new Map();
    const familyMap = new Map();
    const virusMap = new Map();  // Track unique viruses to avoid duplicates
    
    // Count objects to size nodes proportionally
    const phylumCounts = {};
    const familyCounts = {};
    
    // First pass: count items for sizing
    if (data && data.Viral_Hits) {
        data.Viral_Hits.forEach(viralHit => {
            const phylum = viralHit.Phylum || "Unknown Phylum";
            const family = viralHit.Family || "Unknown Family";
            
            // Count phyla and families for sizing nodes
            phylumCounts[phylum] = (phylumCounts[phylum] || 0) + 1;
            
            const familyKey = `${family}_${phylum}`;
            familyCounts[familyKey] = (familyCounts[familyKey] || 0) + 1;
        });
    }
    
    // Second pass: create nodes and links
    if (data && data.Viral_Hits) {
        data.Viral_Hits.forEach(viralHit => {
            const virusName = viralHit.ScientificName || viralHit.QueryID;
            const phylum = viralHit.Phylum || "Unknown Phylum";
            const family = viralHit.Family || "Unknown Family";
            
            // Add phylum if not exists
            if (!phylaMap.has(phylum)) {
                phylaMap.set(phylum, true);
                nodes.push({
                    id: phylum,
                    name: phylum,
                    type: "phylum",
                    count: phylumCounts[phylum] || 1,
                    viralHits: [] // To store refs to viral hits in this phylum
                });
            }
            
            // Add family if not exists
            const familyId = `${family}_${phylum}`;
            if (!familyMap.has(familyId)) {
                familyMap.set(familyId, true);
                nodes.push({
                    id: familyId,
                    name: family,
                    type: "family",
                    phylum: phylum,
                    count: familyCounts[familyId] || 1,
                    viralHits: [] // To store refs to viral hits in this family
                });
                
                // Link family to phylum
                links.push({
                    source: phylum,
                    target: familyId,
                    value: familyCounts[familyId] || 1 // Width based on count
                });
            }
            
            // Add virus (only if unique)
            const virusId = `${virusName}_${familyId}`;
            if (!virusMap.has(virusId)) {
                virusMap.set(virusId, true);
                nodes.push({
                    id: virusId,
                    name: virusName,
                    type: "virus",
                    family: family,
                    phylum: phylum,
                    Genome: viralHit.Genome,
                    taxId: viralHit.TaxId ? Math.trunc(viralHit.TaxId) : null,
                    viralHit: viralHit // Store reference to original data
                });
                
                // Find the family node and add this virus to its viralHits
                const familyNode = nodes.find(node => node.id === familyId);
                if (familyNode) {
                    familyNode.viralHits.push(viralHit);
                }
                
                // Find the phylum node and add this virus to its viralHits
                const phylumNode = nodes.find(node => node.id === phylum);
                if (phylumNode) {
                    phylumNode.viralHits.push(viralHit);
                }
                
                // Link virus to family
                links.push({
                    source: familyId,
                    target: virusId,
                    value: 1
                });
            }
        });
    }
    
    return { nodes, links };
}

// Add search functionality
// Add filter controls instead of just search

}



// call function
populateDashboard(logData);
populateKnownVirusesDashboard(jsonData);
createTaxonomyVisualization(jsonData);


// dealing with JSON dataase
try {
    createIndex(jsonData);
    createVisualizations(jsonData);
} catch (error) {
    console.error("Error loading visualizations:", error);
    document.getElementById("loading").style.display = "none";
    document.getElementById("error").style.display = "block";
}
});


function initializeSelectionAndExport() {
    // CSS styles for new elements
    const styleEl = document.createElement('style');
    styleEl.textContent = `
        .checkbox-container {
            position: absolute; 
            left: 14px;
            top: 14px;
            z-index: 10;
        }
        .sequence-checkbox {
            width: 20px;
            height: 20px;
            cursor: pointer;
        }
        .export-panel {
            margin: 20px 0;
            padding: 15px;
            background-color: #ffffff;
            border: 1px solid #ecf0f1;
            border-radius: 10px;
            display: flex;
            justify-content: space-between;
            box-shadow: 0 2px 4px var(--shadow-color);
        }
        .export-btn {
            background-color: #3498db;
            color: white;
            border: none;
            padding: 10px 15px;
            border-radius: 4px;
            cursor: pointer;
            font-family: system-ui, -apple-system, sans-serif;
            font-size: 16px;
        }
        .export-btn:hover {
            background-color: #2980b9;
        }
        .select-btns {
            display: flex;
            gap: 10px;
        }
        .select-btn {
            background-color: #f1f1f1;
            border: 1px solid #ddd;
            padding: 8px 12px;
            border-radius: 4px;
            cursor: pointer;
            font-family: system-ui, -apple-system, sans-serif;
            font-size: 16px;
        }
        .select-btn:hover {
            background-color: #e4e4e4;
        }
    `;
    document.head.appendChild(styleEl);

    // export panel to the top of visualization container
    const visualizationContainer = document.getElementById('visualization-container');
    const exportPanel = document.createElement('div');
    exportPanel.className = 'export-panel';
    exportPanel.innerHTML = `
        <div class="select-btns">
            <button class="select-btn" id="select-all-btn">Select All</button>
            <button class="select-btn" id="deselect-all-btn">Deselect All</button>
        </div>
        <button class="export-btn" id="export-tsv-btn">Export Selected as TSV</button>
    `;
    visualizationContainer.insertBefore(exportPanel, visualizationContainer.firstChild);

    // checkbox to each visualization wrapper
    const wrappers = document.querySelectorAll('.visualization-wrapper');
    wrappers.forEach(wrapper => {
        const contigId = wrapper.getAttribute('data-contig');
        const checkboxContainer = document.createElement('div');
        checkboxContainer.className = 'checkbox-container';
        
        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.className = 'sequence-checkbox';
        checkbox.setAttribute('data-contig', contigId);
        
        checkboxContainer.appendChild(checkbox);
        wrapper.style.position = 'relative';
        wrapper.appendChild(checkboxContainer);
    });

    // event listeners
    document.getElementById('select-all-btn').addEventListener('click', () => {
        document.querySelectorAll('.sequence-checkbox').forEach(cb => {
            cb.checked = true;
        });
    });

    document.getElementById('deselect-all-btn').addEventListener('click', () => {
        document.querySelectorAll('.sequence-checkbox').forEach(cb => {
            cb.checked = false;
        });
    });

    document.getElementById('export-tsv-btn').addEventListener('click', exportSelectedAsTSV);
}

// export selected data as TSV
function exportSelectedAsTSV() {
    const selected = Array.from(document.querySelectorAll('.sequence-checkbox:checked'))
        .map(cb => cb.getAttribute('data-contig'));
    
    if (selected.length === 0) {
        alert('Please select at least one sequence to export');
        return;
    }
    
    // data for each selected contig
    const exportData = [];
    const headers = [
        'Contig ID', // 'Length', 
        'BLASTx Subject', 'BLASTx Organism', 'BLASTx Coverage', 'BLASTx Identity', 'BLASTx E-value',
        'BLASTn Subject', 'BLASTn Coverage', 'BLASTn Identity', 'BLASTn E-value'
    ];
    
    selected.forEach(contigId => {
        // viral hit data for this contig
        const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${contigId}"]`);
        
        // extract information from the container
        const infoPanel = wrapper.querySelector('.info-panel');
        if (!infoPanel) return;
        
        // parse text content from the info panel
        const infoText = infoPanel.textContent;
        
        // row object
        const row = {
            'Contig ID': contigId
        //'Length': extractValue(infoText, 'Length') || ''
        };
        
        // BLASTx data
        if (infoText.includes('BLASTx Results')) {
            const blastxSection = infoText.split('BLASTx Results')[1].split('BLASTn Results')[0];
            row['BLASTx Subject'] = extractValue(blastxSection, 'Subject') || '';
            row['BLASTx Organism'] = extractValue(blastxSection, 'Organism') || '';
            row['BLASTx Coverage'] = extractValue(blastxSection, 'Coverage') || '';
            row['BLASTx Identity'] = extractValue(blastxSection, 'Identity') || '';
            row['BLASTx E-value'] = extractValue(blastxSection, 'E-value') || '';
        }
        
        // BLASTn data
        if (infoText.includes('BLASTn Results')) {
            const blastnSection = infoText.split('BLASTn Results')[1];
            row['BLASTn Subject'] = extractValue(blastnSection, 'Subject') || '';
            row['BLASTn Coverage'] = extractValue(blastnSection, 'Coverage') || '';
            row['BLASTn Identity'] = extractValue(blastnSection, 'Identity') || '';
            row['BLASTn E-value'] = extractValue(blastnSection, 'E-value') || '';
        }
        
        exportData.push(row);
    });
    
    // TSV content
    let tsvContent = headers.join('\t') + '\n';
    
    exportData.forEach(row => {
        const rowValues = headers.map(header => row[header] || '');
        tsvContent += rowValues.join('\t') + '\n';
    });
    
    // trigger download
    const blob = new Blob([tsvContent], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'viralquest_export.tsv';
    a.click();
    URL.revokeObjectURL(url);
}

// helper function to extract values from text content
function extractValue(text, label) {
    const regex = new RegExp(label + ':\s*([^\n]+)');
    const match = text.match(regex);
    return match ? match[1].trim() : null;
}

// run the initialization after visualizations are created
const observer = new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
        if (mutation.type === 'childList' && 
            document.querySelectorAll('.visualization-wrapper').length > 0 &&
            document.getElementById('loading').style.display === 'none') {
            // visualizations are loaded, initialize our features
            initializeSelectionAndExport();
            // stop observing once initialization is complete
            observer.disconnect();
        }
    });
});

// start observing
observer.observe(document.getElementById('visualization-container'), { childList: true, subtree: true });

// also try to initialize after a few seconds
setTimeout(function() {
    if (document.querySelectorAll('.visualization-wrapper').length > 0) {
        initializeSelectionAndExport();
        observer.disconnect();
    }
}, 3000);

// initialize enhanced export functionality
function initializeEnhancedExport() {
// export panel HTML to include format options
const exportPanel = document.querySelector('.export-panel');
exportPanel.innerHTML = `
    <div class="select-btns">
    <button class="select-btn" id="select-all-btn">Select All</button>
    <button class="select-btn" id="deselect-all-btn">Deselect All</button>
    </div>
    <div class="export-options">
    <button class="export-btn" id="export-tsv-btn">Export Selected</button>
    <div class="export-format">
        <label>Format:</label>
        <select id="export-format">
        <option value="tsv">TSV (Data)</option>
        <option value="svg">SVG (Graphics)</option>
        </select>
    </div>
    <div class="toggle-columns" id="customize-columns">Customize Columns</div>
    </div>
`;

// columns dialog
const columnsDialog = document.createElement('div');
columnsDialog.className = 'columns-dialog';
columnsDialog.id = 'columns-dialog';

// overlay
const overlay = document.createElement('div');
overlay.className = 'overlay';
overlay.id = 'dialog-overlay';

document.body.appendChild(columnsDialog);
document.body.appendChild(overlay);

// event listeners
document.getElementById('select-all-btn').addEventListener('click', () => {
    document.querySelectorAll('.sequence-checkbox').forEach(cb => {
    cb.checked = true;
    });
});

document.getElementById('deselect-all-btn').addEventListener('click', () => {
    document.querySelectorAll('.sequence-checkbox').forEach(cb => {
    cb.checked = false;
    });
});

document.getElementById('export-tsv-btn').addEventListener('click', () => {
    const format = document.getElementById('export-format').value;
    if (format === 'tsv') {
    exportSelectedAsTSV();
    } else if (format === 'svg') {
    exportSelectedAsGraphics('svg');
    } else if (format === 'pdf') {
    exportSelectedAsGraphics('pdf');
    }
});

// customization dialog
document.getElementById('customize-columns').addEventListener('click', showColumnsDialog);

// close dialog when clicking overlay
document.getElementById('dialog-overlay').addEventListener('click', hideColumnsDialog);
}

// columns selection dialog
function showColumnsDialog() {
    const columnsDialog = document.getElementById('columns-dialog');
    const overlay = document.getElementById('dialog-overlay');

// define all possible columns from jsonDB
const possibleColumns = [
    // viral Hits columns
    { id: 'Sample_name', label: 'Sample Name', category: 'General', checked: true },
    { id: 'QueryID', label: 'Contig ID', category: 'General', checked: true },
    { id: 'BLASTx_Qlength', label: 'BLASTx Query Length', category: 'BLASTx', checked: true },
    { id: 'BLASTx_Slength', label: 'BLASTx Subj. Length', category: 'BLASTx', checked: true },
    { id: 'BLASTx_Cover', label: 'BLASTx Coverage', category: 'BLASTx', checked: true },
    { id: 'BLASTx_Ident', label: 'BLASTx Identity', category: 'BLASTx', checked: true },
    { id: 'BLASTx_evalue', label: 'BLASTx E-value', category: 'BLASTx', checked: true },
    { id: 'BLASTx_Subject_Title', label: 'BLASTx Subject Title', category: 'BLASTx', checked: true },
    { id: 'BLASTx_Organism_Name', label: 'BLASTx Organism', category: 'BLASTx', checked: true },
    { id: 'BLASTn_Qlength', label: 'BLASTn Query Length', category: 'BLASTn', checked: false },
    { id: 'BLASTn_Slength', label: 'BLASTn Subj. Length', category: 'BLASTn', checked: false },
    { id: 'BLASTn_Cover', label: 'BLASTn Coverage', category: 'BLASTn', checked: false },
    { id: 'BLASTn_Ident', label: 'BLASTn Identity', category: 'BLASTn', checked: false },
    { id: 'BLASTn_evalue', label: 'BLASTn E-value', category: 'BLASTn', checked: false },
    { id: 'BLASTn_Subject_Title', label: 'BLASTn Subject Title', category: 'BLASTn', checked: false },
    { id: 'TaxId', label: 'Taxonomy ID', category: 'Taxonomy', checked: false },
    { id: 'ScientificName', label: 'Scientific Name', category: 'Taxonomy', checked: true },
    { id: 'Clade', label: 'Clade', category: 'Taxonomy', checked: false },
    { id: 'Kingdom', label: 'Kingdom', category: 'Taxonomy', checked: true },
    { id: 'Phylum', label: 'Phylum', category: 'Taxonomy', checked: true },
    { id: 'Class', label: 'Class', category: 'Taxonomy', checked: true },
    { id: 'Order', label: 'Order', category: 'Taxonomy', checked: true },
    { id: 'Family', label: 'Family', category: 'Taxonomy', checked: true },
    { id: 'Genus', label: 'Genus', category: 'Taxonomy', checked: true },
    { id: 'Species', label: 'Species', category: 'Taxonomy', checked: true },
    { id: 'Genome', label: 'Genome Type', category: 'Taxonomy', checked: true }
];

// group columns by category
const categories = [...new Set(possibleColumns.map(col => col.category))];

// dialog content
let dialogContent = `
    <h3>Customize Export Columns</h3>
    <p>Select the columns to include in your TSV export:</p>
`;

categories.forEach(category => {
    dialogContent += `<h4>${category}</h4>`;
    dialogContent += '<div class="columns-grid">';
    
    possibleColumns.filter(col => col.category === category).forEach(col => {
        dialogContent += `
            <div class="column-item">
            <input type="checkbox" id="col-${col.id}" data-column="${col.id}" ${col.checked ? 'checked' : ''}>
            <label for="col-${col.id}">${col.label}</label>
            </div>
        `;
    });
    
    dialogContent += '</div>';
});

dialogContent += `
    <div class="columns-dialog-buttons">
    <button id="select-all-columns">Select All</button>
    <button id="deselect-all-columns">Deselect All</button>
    <button id="apply-columns">Apply</button>
    <button id="cancel-columns">Cancel</button>
    </div>
`;

columnsDialog.innerHTML = dialogContent;

// event listeners to dialog buttons
document.getElementById('select-all-columns').addEventListener('click', () => {
    columnsDialog.querySelectorAll('input[type="checkbox"]').forEach(cb => {
        cb.checked = true;
    });
});    

document.getElementById('deselect-all-columns').addEventListener('click', () => {
    columnsDialog.querySelectorAll('input[type="checkbox"]').forEach(cb => {
        cb.checked = false;
    });
});

document.getElementById('apply-columns').addEventListener('click', () => {
    // save selected columns to localStorage for persistence
    const selectedColumns = Array.from(columnsDialog.querySelectorAll('input[type="checkbox"]'))
        .filter(cb => cb.checked)
        .map(cb => cb.dataset.column);
    
    localStorage.setItem('vqSelectedColumns', JSON.stringify(selectedColumns));
    hideColumnsDialog();
});

document.getElementById('cancel-columns').addEventListener('click', hideColumnsDialog);

// dialog and overlay
columnsDialog.style.display = 'block';
overlay.style.display = 'block';
}

// hide columns selection dialog
function hideColumnsDialog() {
    document.getElementById('columns-dialog').style.display = 'none';
    document.getElementById('dialog-overlay').style.display = 'none';
}

// function to export selected contigs as TSV with data from jsonDB
function exportSelectedAsTSV() {
    const selected = Array.from(document.querySelectorAll('.sequence-checkbox:checked'))
        .map(cb => cb.getAttribute('data-contig'));
    
    if (selected.length === 0) {
        alert('Please select at least one sequence to export');
        return;
    }
    
    // selected columns from localStorage or use defaults
    let selectedColumns;
    try {
        selectedColumns = JSON.parse(localStorage.getItem('vqSelectedColumns')) || [
        'QueryID', 'Sample_name', 'BLASTx_Qlength', 'BLASTx_Cover', 'BLASTx_Ident', 
        'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name',
        'ScientificName', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
        ];
    } catch (e) {
        // default columns if parsing fails
        selectedColumns = [
        'QueryID', 'Sample_name', 'BLASTx_Qlength', 'BLASTx_Cover', 'BLASTx_Ident', 
        'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name',
        'ScientificName', 'Kingdom', 'Family', 'Genus', 'Species'
        ];
    }
    
    // column labels (for header row)
    const columnLabels = {
        'QueryID': 'Contig ID',
        'Sample_name': 'Sample_Name',
        'BLASTx_Qlength': 'BLASTx_Length',
        'BLASTx_Slength': 'BLASTx_Subject_Length',
        'BLASTx_Cover': 'BLASTx_Coverage',
        'BLASTx_Ident': 'BLASTx_Identity',
        'BLASTx_evalue': 'BLASTx_Evalue',
        'BLASTx_Subject_Title': 'BLASTx_Subject',
        'BLASTx_Organism_Name': 'BLASTx_Organism',
        'BLASTn_Qlength': 'BLASTn_Length',
        'BLASTn_Slength': 'BLASTn_Subject_Length',
        'BLASTn_Cover': 'BLASTn_Coverage',
        'BLASTn_Ident': 'BLASTn_Identity',
        'BLASTn_evalue': 'BLASTn_Evalue',
        'BLASTn_Subject_Title': 'BLASTn_Subject',
        'TaxId': 'Taxonomy_ID',
        'ScientificName': 'Scientific_Name',
        'Clade': 'Clade',
        'Kingdom': 'Kingdom',
        'Phylum': 'Phylum',
        'Class': 'Class',
        'Order': 'Order',
        'Family': 'Family',
        'Subfamily': 'Subfamily',
        'Genus': 'Genus',
        'Species': 'Species',
        'Genome': 'Genome Type'
    };
    
    // header row
    const headers = selectedColumns.map(col => columnLabels[col] || col);
    
    // data for selected contigs directly from jsonDB
    const rows = [];
    
    selected.forEach(contigId => {
        // find contigs Viral_Hits data
        const hit = jsonData.Viral_Hits.find(h => h.QueryID === contigId);
        if (hit) {
        const row = {};
        
        // values for selected columns
        selectedColumns.forEach(col => {
            row[columnLabels[col] || col] = hit[col] !== undefined ? hit[col] : '';
        });
        
        rows.push(row);
        }
    });
    
    // TSV content
    let tsvContent = headers.join('\t') + '\n';
    
    rows.forEach(row => {
        const rowValues = headers.map(header => row[header] || '');
        tsvContent += rowValues.join('\t') + '\n';
    });
    
    // trigger download
    const blob = new Blob([tsvContent], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'viralquest_export.tsv';
    a.click();
    URL.revokeObjectURL(url);
}

// export selected visualizations as SVG or PDF
function exportSelectedAsGraphics(format) {
    const selected = Array.from(document.querySelectorAll('.sequence-checkbox:checked'))
        .map(cb => cb.getAttribute('data-contig'));
    
    if (selected.length === 0) {
        alert('Please select at least one visualization to export');
        return;
    }
    
    if (format === 'svg') {
        // Export as individual SVG files or as a combined SVG
        if (selected.length === 1) {
        // Single SVG export
        exportSingleSVG(selected[0]);
        } else {
        // multiple SVGs - ask if user wants individual files or a combined file
        if (confirm('Export as individual SVG files? Click Cancel for a combined file.')) {
            selected.forEach(contigId => exportSingleSVG(contigId));
        } else {
            exportCombinedSVG(selected);
        }
        }
    } else if (format === 'pdf') {
        // for PDF we need to use a library like jsPDF
        // check if jsPDF is loaded
        if (typeof jspdf === 'undefined') {
        // load jsPDF dynamically
        const script = document.createElement('script');
        script.src = 'https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js';
        script.onload = function() {
            // also load svg2pdf
            const svg2pdf = document.createElement('script');
            svg2pdf.src = 'https://cdnjs.cloudflare.com/ajax/libs/svg2pdf.js/2.2.1/svg2pdf.min.js';
            svg2pdf.onload = function() {
            exportAsPDF(selected);
            };
            document.head.appendChild(svg2pdf);
        };
        document.head.appendChild(script);
        } else {
        // jsPDF is already loaded
        exportAsPDF(selected);
        }
    }
}

// wxport a single visualization as SVG
function exportSingleSVG(contigId) {
    const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${contigId}"]`);
    if (!wrapper) return;
    
    const svg = wrapper.querySelector('svg');
    if (!svg) return;
    
    // clone the SVG to make modifications without affecting the original
    const svgClone = svg.cloneNode(true);
    
    // optional: clean up or modify SVG for export
    enhanceSVGForExport(svgClone);
    
    // SVG content
    const svgData = new XMLSerializer().serializeToString(svgClone);
    
    // create a downloadable blob
    const blob = new Blob([svgData], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    
    // trigger download
    const a = document.createElement('a');
    a.href = url;
    a.download = `${contigId}.svg`;
    a.click();
    URL.revokeObjectURL(url);
}

// export multiple visualizations as a combined SVG
function exportCombinedSVG(contigIds) {
    // new SVG to hold all the visualizations
    const combinedSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    
    // initialize height and width tracking variables
    let totalHeight = 20; // Start with margin
    let maxWidth = 0;
    
    // container for the title
    const titleG = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    
    // title text
    const titleText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    titleText.setAttribute('x', 20);
    titleText.setAttribute('y', 20);
    titleText.setAttribute('font-family', 'Arial, sans-serif');
    titleText.setAttribute('font-size', '16px');
    titleText.setAttribute('font-weight', 'bold');
    //titleText.textContent = 'ViralQuest Visualization Export';
    
    titleG.appendChild(titleText);
    combinedSvg.appendChild(titleG);
    
    totalHeight += 30; // Add space after title
    
    // process each selected contig
    contigIds.forEach((contigId, index) => {
        const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${contigId}"]`);
        if (!wrapper) return;
        
        const svg = wrapper.querySelector('svg');
        if (!svg) return;
        
        // visualization content 
        const mainGroup = svg.querySelector('g');
        if (!mainGroup) return;
        
        // clone it
        const groupClone = mainGroup.cloneNode(true);
        
        // position the group at the correct vertical position
        groupClone.setAttribute('transform', `translate(20, ${totalHeight})`);
        
        // label for this contig
        const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        label.setAttribute('x', 0);
        label.setAttribute('y', -10);
        label.setAttribute('font-family', 'Arial, sans-serif');
        label.setAttribute('font-size', '14px');
        label.setAttribute('font-weight', 'bold');
        label.textContent = contigId;
        
        groupClone.insertBefore(label, groupClone.firstChild);
        
        // combined SVG
        combinedSvg.appendChild(groupClone);
        
        // update height and width tracking
        const bbox = mainGroup.getBBox();
        totalHeight += bbox.height + 50; // Add space between contigs
        maxWidth = Math.max(maxWidth, bbox.width + 40); // Add margin
    });
    
    // set dimensions of the combined SVG
    combinedSvg.setAttribute('width', maxWidth);
    combinedSvg.setAttribute('height', totalHeight);
    combinedSvg.setAttribute('viewBox', `0 0 ${maxWidth} ${totalHeight}`);
    
    // styling
    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
    style.textContent = `
        text { font-family: system-ui, -apple-system, sans-serif; }
        .genome-line { stroke: #2c3e50; stroke-width: 2; }
        .scale-axis path, .scale-axis line { stroke: #7f8c8d; }
    `;
    combinedSvg.insertBefore(style, combinedSvg.firstChild);
    
    // SVG content
    const svgData = new XMLSerializer().serializeToString(combinedSvg);
    
    // downloadable blob
    const blob = new Blob([svgData], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    
    // trigger download
    const a = document.createElement('a');
    a.href = url;
    a.download = 'viralquest_combined.svg';
    a.click();
    URL.revokeObjectURL(url);
}

// enhance SVG before export
function enhanceSVGForExport(svgElement) {
    // CSS styles directly to the SVG
    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
    style.textContent = `
        text { font-family: system-ui, -apple-system, sans-serif; }
        .genome-line { stroke: #2c3e50; stroke-width: 2; }
        .scale-axis path, .scale-axis line { stroke: #7f8c8d; }
        .pfam-domain { stroke: #333; stroke-width: 0.5; }
    `;
    svgElement.insertBefore(style, svgElement.firstChild);
    
    // title and description
    const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
    title.textContent = 'ViralQuest Visualization';
    svgElement.insertBefore(title, svgElement.firstChild);
    
    // ensure elements have proper styling attributes
    const paths = svgElement.querySelectorAll('path');
    paths.forEach(path => {
        if (!path.getAttribute('stroke') && !path.getAttribute('style').includes('stroke')) {
        path.setAttribute('stroke', '#333');
        path.setAttribute('stroke-width', '0.5');
        }
    });
    
    return svgElement;
}




function enhanceExistingExport() {
    // check if basic export panel exists
    const existingPanel = document.querySelector('.export-panel');
    if (existingPanel) {
        // replace with enhanced version
        initializeEnhancedExport();
    } else {
        setTimeout(enhanceExistingExport, 500);
    }
}

// call function
document.addEventListener('DOMContentLoaded', () => {
    // check basic export panel
    setTimeout(enhanceExistingExport, 1000);
});
