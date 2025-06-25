import os, json

from .info import __version__

def merge_json_files(directory_path: str) -> None:
    """
    merge multiple JSON files with the same structure into a single JSON file.
    """
    if not os.path.isdir(directory_path):
        raise ValueError(f"The provided path '{directory_path}' is not a valid directory")
    
    output_file = "ViralQuestReport.html"
    
    # Initialize the merged data structure
    merged_data = {
        "Viral_Hits": [],
        "HMM_hits": [],
        "ORF_Data": {}
    }
    
    # Get all JSON files in the directory
    json_files = [f for f in os.listdir(directory_path) if f.endswith('.json')]
    
    if not json_files:
        print(f"No JSON files found in {directory_path}")
        return
    
    print(f"Found {len(json_files)} JSON files to merge")
    
    # Process each JSON file
    for file_name in json_files:
        file_path = os.path.join(directory_path, file_name)
        try:
            with open(file_path, 'r') as file:
                data = json.load(file)
                
                # Merge Viral_Hits (list of dictionaries)
                if "Viral_Hits" in data:
                    merged_data["Viral_Hits"].extend(data["Viral_Hits"])
                
                # Merge HMM_hits (list of dictionaries)
                if "HMM_hits" in data:
                    merged_data["HMM_hits"].extend(data["HMM_hits"])
                
                # Merge ORF_Data (dictionary of lists)
                if "ORF_Data" in data:
                    for contig_id, orf_list in data["ORF_Data"].items():
                        if contig_id in merged_data["ORF_Data"]:
                            # If the contig ID already exists, extend its list
                            merged_data["ORF_Data"][contig_id].extend(orf_list)
                        else:
                            # If the contig ID is new, add it to the merged data
                            merged_data["ORF_Data"][contig_id] = orf_list
                
            #print(f"Successfully merged data from {file_name}")
        except json.JSONDecodeError:
            print(f"Error: {file_name} is not a valid JSON file. Skipping.")
        except Exception as e:
            print(f"Error processing {file_name}: {str(e)}. Skipping.")
    
        jsonDB = merged_data
        with open(output_file, 'w') as outfile:
            json.dump(jsonDB, outfile, indent=4)

    html_content = fr"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>ViralQuest</title>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
        <style>
            :root {{
                --primary-color: #2c3e50;
                --secondary-color: #3498db; 
                --background-color: #f5f7fa;
                --card-background: #ffffff;
                --border-color: #e2e8f0;
                --shadow-color: rgba(0, 0, 0, 0.1);
                --text-color: #2d3748;
                --accent-color: #4299e1;
                --pfam-color-1: #ff6b6b;
                --pfam-color-2: #58b368;
                --pfam-color-3: #bc6ff1;
                --pfam-color-4: #ffa048;
                --pfam-color-5: #5191d1;
                --pfam-color-6: #ffcc29;
                --pfam-color-7: #f06595;
                --pfam-color-8: #38b2ac;
                --success-color: #2ecc71;
                --warning-color: #f39c12;
                --card-background: #ffffff;
                --border-color: #e0e0e0;
                --text-primary: #2c3e50;
                --text-secondary: #7f8c8d;
                --chart-color-1: #3498db;
                --chart-color-2: #2ecc71;
                --chart-color-3: #e74c3c;
            }}

            body {{
                margin: 0;
                padding: 0;
                background-color: var(--background-color);
                color: var(--text-color);
                font-family: system-ui, -apple-system, sans-serif;
            }}

            /* container styles */
            .page-container {{
                max-width: 1400px;
                margin: 0 auto;
                padding: 20px;
                display: flex;
                gap: 20px;
                position: relative;
            }}


            header {{
                background: var(--card-background);
                padding: 20px;
                border-radius: 12px;
                box-shadow: 0 2px 4px var(--shadow-color);
                margin-bottom: 20px;
            }}

            header h1 {{
                color: var(--primary-color);
                margin: 0 0 10px 0;
                font-size: 2.5rem;
            }}

            header h3 {{
                color: var(--text-color);
                margin: 0;
                font-weight: 500;
                opacity: 0.8;
            }}

            /* index styles */
            #contig-index {{
                width: 10%;
                background: var(--card-background);
                border-radius: 12px;
                padding: 20px;
                box-shadow: 0 2px 4px var(--shadow-color);
                position: fixed;
                top: 1%;  
                bottom: 1%;
                left: 20px;
                /*max-height: calc(100vh - 140px);*/
                display: flex;
                flex-direction: column;
                transition: transform 0.3s ease;
                z-index: 1000;
            }}



            #contig-index h2 {{
                color: var(--primary-color);
                margin: 0 0 15px 0;
                font-size: 1.2rem;
                padding-right: 30px; 
            }}

            .index-list {{
                list-style: none;
                padding: 0;
                margin: 0;
                overflow-y: auto;
                flex-grow: 1;
                scrollbar-width: thin;
                scrollbar-color: var(--primary-color) var(--background-color);
            }}

            .index-item {{
                padding: 8px 12px;
                margin: 4px 0;
                border-radius: 6px;
                cursor: pointer;
                transition: all 0.2s ease;
                color: var(--text-color);
                font-size: 0.9rem;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
            }}


            .page-container {{
                margin-left: 290px; /* width of index + some margin */
                max-width: calc(100% - 310px); /* adjust based on index width */
                transition: margin-left 0.3s ease;
                margin-right: 20%;  
            }}

            .index-item:hover {{
                background: var(--accent-color);
                color: white;
            }}

            /* responsive adjustments */
            @media (max-width: 768px) {{
                #contig-index {{
                    width: 200px;
                }}

                .page-container {{
                    margin-left: 240px;
                    max-width: calc(100% - 260px);
                }}
            }}

            /* main visualization container */
            #visualization-container {{
                flex-grow: 1;
                background: transparent;
            }}

            /* individual visualization wrapper */
            .visualization-wrapper {{
                background: var(--card-background);
                margin-bottom: 30px;
                border-radius: 12px;
                padding: 25px;
                box-shadow: 0 4px 6px var(--shadow-color);
                transition: transform 0.2s ease;
            }}

            .visualization-wrapper:hover {{
                transform: translateY(-2px);
            }}

            /* tooltip styles */
            [class^='tooltip-'] {{
                position: absolute;
                display: none;
                background: var(--card-background);
                border: 1px solid var(--border-color);
                border-radius: 8px;
                padding: 15px;
                box-shadow: 0 4px 12px var(--shadow-color);
                max-width: 350px;
                z-index: 1000;
                font-size: 0.9rem;
            }}

            /* Copy button */
            .copy-button {{
                margin-left: 8px;
                padding: 6px 12px;
                font-size: 0.8rem;
                cursor: pointer;
                border: 1px solid var(--border-color);
                background: var(--card-background);
                border-radius: 4px;
                transition: all 0.2s ease;
                color: var(--primary-color);
            }}

            .copy-button:hover {{
                background: var(--accent-color);
                color: white;
                border-color: var(--accent-color);
            }}

            .copy-feedback {{
                display: none;
                margin-left: 8px;
                color: #10b981;
                font-size: 0.8rem;
                opacity: 0;
                transition: opacity 0.3s ease;
            }}

            .copy-feedback.show {{
                opacity: 1;
            }}

            /* Loading and error */
            .loading, .error {{
                text-align: center;
                padding: 30px;
                font-size: 1rem;
                border-radius: 8px;
                margin: 20px 0;
            }}

            .loading {{
                background: var(--card-background);
                color: var(--text-color);
                box-shadow: 0 2px 4px var(--shadow-color);
            }}

            .error {{
                background: #fee2e2;
                color: #dc2626;
                border: 1px solid #fecaca;
            }}

            /* Pfam domain */
            .pfam-domain {{
                cursor: pointer;
                stroke: #ffffff;
                stroke-width: 1px;
                transition: all 0.2s ease;
            }}

            .pfam-domain:hover {{
                stroke-width: 2px;
                filter: brightness(1.1);
            }}

            .pfam-tooltip {{
                max-width: 400px;
                max-height: 300px;
                overflow-y: auto;
            }}

            .pfam-info {{
                font-size: 0.85rem;
                line-height: 1.4;
                margin-top: 8px;
                padding: 8px;
                border-radius: 4px;
                background-color: #f8fafc;
                border-left: 3px solid var(--accent-color);
            }}

            /* scrollbar styling */
            ::-webkit-scrollbar {{
                width: 8px;
            }}

            ::-webkit-scrollbar-track {{
                background: var(--background-color);
                border-radius: 4px;
            }}

            ::-webkit-scrollbar-thumb {{
                background: var(--primary-color);
                border-radius: 4px;
            }}

            ::-webkit-scrollbar-thumb:hover {{
                background: var(--accent-color);
            }}

            .stats-dashboard {{
                background-color: #f9f9f9;
                border-radius: 8px;
                margin: 1% auto;
                margin-right: 15.5%;
                padding: 1%;
                max-width: 65%;
                box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
            }}

            .stats-dashboard h2 {{
                color: var(--text-primary);
                margin-top: 0;
                margin-bottom: 20px;
                font-size: 24px;
                text-align: center;
            }}

            .dashboard-container {{
                display: grid;
                grid-template-columns: 1fr;
                gap: 20px;
            }}

            @media (min-width: 768px) {{
                .dashboard-container {{
                    grid-template-columns: 1fr 1fr;
                }}
            }}

            .stats-cards {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
            }}

            .stat-card {{
                background-color: var(--card-background);
                border-radius: 8px;
                padding: 15px;
                display: flex;
                align-items: center;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.05);
                transition: transform 0.2s ease, box-shadow 0.2s ease;
            }}

            .stat-card:hover {{
                transform: translateY(-3px);
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            }}

            .stat-card.highlight {{
                border-left: 4px solid var(--accent-color);
            }}

            .stat-icon {{
                display: flex;
                align-items: center;
                justify-content: center;
                width: 48px;
                height: 48px;
                border-radius: 50%;
                background-color: rgba(52, 152, 219, 0.1);
                color: var(--primary-color);
                margin-right: 15px;
            }}

            .stat-content {{
                flex: 1;
            }}

            .stat-value {{
                font-size: 24px;
                font-weight: bold;
                color: var(--text-primary);
                margin-bottom: 5px;
            }}

            .stat-label {{
                font-size: 14px;
                color: var(--text-secondary);
            }}

            .charts-container {{
                display: grid;
                grid-template-columns: 1fr;
                gap: 20px;
            }}

            .chart-box {{
                background-color: var(--card-background);
                border-radius: 8px;
                padding: 20px;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.05);
            }}

            .chart-box h3 {{
                margin-top: 0;
                margin-bottom: 15px;
                color: var(--text-primary);
                font-size: 18px;
            }}

            #sequences-chart {{
                height: 200px;
                margin-top: 10px;
            }}

            .config-info {{
                display: flex;
                flex-direction: column;
                gap: 10px;
            }}

            .config-item {{
                display: flex;
                flex-direction: column;
            }}

            .config-label {{
                font-weight: bold;
                color: var(--text-primary);
                font-size: 14px;
                margin-bottom: 3px;
            }}

            .config-value {{
                color: var(--text-secondary);
                word-break: break-word;
                font-family: monospace;
                font-size: 14px;
            }}

            .db-list {{
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                margin-top: 5px;
            }}

            .db-item {{
                background-color: rgba(52, 152, 219, 0.1);
                color: var(--primary-color);
                padding: 5px 10px;
                border-radius: 4px;
                font-size: 13px;
            }}

            /* Animation for loading elements */
            @keyframes fadeIn {{
                from {{ opacity: 0; }}
                to {{ opacity: 1; }}
            }}

            .stat-value, .config-value, .db-item {{
                animation: fadeIn 0.5s ease-in-out;
            }}

            .version_text {{
                font-size: 12px;
                text-align: center;
            }}
            
            .index-logo {{
                width: 100%;
                height: auto;
                display: block;
                margin-bottom: 0px;
                pointer-events: auto;
                cursor: pointer;
            }}

            #contig-index a {{
                position: relative;
                z-index: 1001;
                display: inline-block;
            }}

            .logo-link {{
                display: block;
                margin-bottom: 0px;
                cursor: pointer;
                z-index: 1001;
            }}

            .export-options {{
              display: flex;
              gap: 10px;
              align-items: center;
            }}

            .export-format {{
              display: flex;
              align-items: center;
              margin-left: 15px;
              border-left: 1px solid #ddd;
              padding-left: 15px;
            }}

            .export-format select {{
              margin-left: 8px;
              padding: 8px;
              border-radius: 4px;
              border: 1px solid #ddd;
              font-family: system-ui, -apple-system, sans-serif;
            }}

            .toggle-columns {{
              margin-left: 15px;
              cursor: pointer;
              color: #3498db;
              font-size: 14px;
              text-decoration: underline;
            }}

            .columns-dialog {{
              position: fixed;
              top: 50%;
              left: 50%;
              transform: translate(-50%, -50%);
              background-color: white;
              padding: 20px;
              border-radius: 8px;
              box-shadow: 0 0 20px rgba(0,0,0,0.3);
              z-index: 1000;
              max-width: 600px;
              max-height: 80vh;
              overflow-y: auto;
              display: none;
            }}

            .columns-grid {{
              display: grid;
              grid-template-columns: repeat(3, 1fr);
              gap: 10px;
              margin-top: 15px;
            }}

            .column-item {{
              display: flex;
              align-items: center;
            }}

            .columns-dialog-buttons {{
              margin-top: 20px;
              display: flex;
              justify-content: flex-end;
              gap: 10px;
            }}

            .columns-dialog-buttons button {{
              padding: 8px 15px;
              border-radius: 4px;
              cursor: pointer;
              font-family: system-ui, -apple-system, sans-serif;
            }}

            .overlay {{
              position: fixed;
              top: 0;
              left: 0;
              right: 0;
              bottom: 0;
              background-color: rgba(0,0,0,0.5);
              z-index: 999;
              display: none;
            }}

            .known-viruses-dashboard {{
                background-color: #f9f9f9;
                border-radius: 8px;
                margin: 1% auto;
                margin-right: 15.5%;
                padding: 1%;
                max-width: 65%;
                box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
            }}

            .known-viruses-description {{
                text-align: center;
                margin-bottom: 20px;
                color: var(--text-secondary);
                font-size: 14px;
            }}

            .known-viruses-container {{
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
                gap: 20px;
            }}

            .known-virus-card {{
                background-color: var(--card-background);
                border-radius: 8px;
                padding: 20px;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.05);
                transition: transform 0.2s ease, box-shadow 0.2s ease;
                border-left: 4px solid var(--accent-color);
                display: flex;
                flex-direction: column;
            }}

            .known-virus-card:hover {{
                transform: translateY(-3px);
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            }}

            .virus-header {{
                display: flex;
                align-items: center;
                margin-bottom: 15px;
            }}

            .virus-icon {{
                display: flex;
                align-items: center;
                justify-content: center;
                width: 36px;
                height: 36px;
                border-radius: 50%;
                background-color: rgba(52, 152, 219, 0.1);
                color: var(--primary-color);
                margin-right: 15px;
            }}

            .virus-name {{
                font-size: 16px;
                font-weight: bold;
                color: var(--text-primary);
                flex: 1;
                word-break: break-word;
            }}

            .virus-metrics {{
                display: flex;
                margin-top: 10px;
                justify-content: space-between;
            }}

            .virus-metric {{
                text-align: center;
                flex: 1;
            }}

            .metric-value {{
                font-size: 18px;
                font-weight: bold;
                color: var(--text-primary);
            }}

            .metric-label {{
                font-size: 12px;
                color: var(--text-secondary);
            }}

            .virus-taxonomy {{
                margin-top: 15px;
                padding-top: 10px;
                border-top: 1px solid #ecf0f1;
                font-size: 13px;
                color: var(--text-secondary);
            }}

            .virus-taxonomy span {{
                display: inline-block;
                margin-right: 8px;
                color: var(--text-primary);
            }}

            .contig-link {{
                margin-top: 15px;
                font-size: 13px;
                color: var(--primary-color);
                text-decoration: none;
                align-self: flex-end;
                cursor: pointer;
            }}

            .contig-link:hover {{
                text-decoration: underline;
            }}

            .no-known-viruses {{
                grid-column: 1 / -1;
                text-align: center;
                color: var(--text-secondary);
                padding: 30px;
                font-style: italic;
            }}

            .taxonomy-container {{
                grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
                gap: 20px;
            }}

        </style>
    </head>
    <body>

        <!-- Taxonomy graph -->
        <div class="stats-dashboard known-viruses-dashboard">
            <h2>Viral Taxonomy</h2>
            <div class="known-viruses-description">
                <p>Network visualization of viruses grouped by phylum and family</p>
            </div>
            <div class="taxonomy-container" id="taxonomy-container">
                <!-- Known virus cards will be dynamically added here -->
                <div class="no-known-viruses" id="no-known-viruses">
                    No viruses meet the criteria (≥90% identity & ≥70% coverage)
                </div>
            </div>
        </div>

        <!-- Known Viruses Dashboard to be placed after the existing stats-dashboard div -->
        <div class="stats-dashboard known-viruses-dashboard">
            <h2>Known Viruses</h2>
            <div class="known-viruses-description">
                <p>Identified viral sequences with high confidence (≥90% identity & ≥70% coverage in BLASTx | ≥90% identity/coverage in BLASTn)</p>
            </div>
            <div class="known-viruses-container" id="known-viruses-container">
                <!-- Known virus cards will be dynamically added here -->
                <div class="no-known-viruses" id="no-known-viruses">
                    No viruses meet the criteria (≥90% identity & ≥70% coverage)
                </div>
            </div>
        </div>

        <div class="page-container">
            <aside id="contig-index">
                <a href="https://github.com/gabrielvpina/viralquest" target="_blank" class="logo-link">
                	<img class="index-logo" src="data:image/png;base64,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" alt="Logo">
                </a>
                <p class="version_text"> Version {__version__} </p>
                <h2>Sequences Index</h2>
                <ul class="index-list" id="index-list">
                    <!-- Index items will be dynamically added here -->
                </ul>
            </aside>

            <div id="visualization-container">
                <div id="loading" class="loading">Loading visualizations...</div>
                <div id="error" class="error" style="display: none;">
                    Error loading visualizations. Please try again.
                </div>
            </div>
        </div>

        <script>

        const jsonData = {jsonDB}

        function populateKnownVirusesDashboard(data) {{
            const container = document.getElementById('known-viruses-container');
            const noKnownVirusesElement = document.getElementById('no-known-viruses');
            
            // Filter viral hits that meet the criteria (≥90% identity, ≥70% coverage in BLASTx)
            const knownViruses = data.Viral_Hits.filter(hit => {{
            // condition 1: BLASTx criteria
            const isKnownByBlastX = hit.BLASTx_Ident >= 90 && hit.BLASTx_Cover >= 70;

            // condition 2: BLASTn criteria
            const subjectTitle = (hit.BLASTn_Subject_Title || '').toLowerCase();
            const hasKeyword = subjectTitle.includes('virus') ||
                                subjectTitle.includes('phage') ||
                                subjectTitle.includes('riboviria');
            const isKnownByBlastN = hasKeyword && hit.BLASTn_Ident >= 90 && hit.BLASTn_Cover >= 90;

            // return true if EITHER condition is met
            return isKnownByBlastX || isKnownByBlastN;
            }});
            
            // Show or hide the "no known viruses" message
            if (knownViruses.length === 0) {{
                noKnownVirusesElement.style.display = 'block';
                return;
            }} else {{
                noKnownVirusesElement.style.display = 'none';
            }}
            
            // Create a card for each known virus
            knownViruses.forEach(virus => {{
                const card = document.createElement('div');
                card.className = 'known-virus-card';
                
                // Virus icon (you can use the same svg as in your existing cards)
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
                
                // Create the card content
                card.innerHTML = `
                    <div class="virus-header">
                        <div class="virus-icon">${{virusIcon}}</div>
                        <div class="virus-name">${{organismName}}</div>
                    </div>
                    
                    <div class="virus-metrics">
                        <div class="virus-metric">
                            <div class="metric-value">${{virus.BLASTx_Ident}}%</div>
                            <div class="metric-label">Identity</div>
                        </div>
                        <div class="virus-metric">
                            <div class="metric-value">${{virus.BLASTx_Cover}}%</div>
                            <div class="metric-label">Coverage</div>
                        </div>
                        <div class="virus-metric">
                            <div class="metric-value">${{Math.trunc(virus.BLASTx_Qlength)}} nt</div>
                            <div class="metric-label">Length</div>
                        </div>
                    </div>
                    
                    <div class="virus-taxonomy">
                        ${{virus.Family ? `<span>Family: ${{virus.Family}}</span>` : ''}}
                        ${{virus.Genome ? `<span>Genome: ${{virus.Genome}}</span>` : ''}}
                    </div>
                    
                    <a class="contig-link" onclick="scrollToContig('${{virus.QueryID}}')">
                        View details (${{virus.QueryID}})
                    </a>
                `;
                
                container.appendChild(card);
            }});
            
            // Update the stats in the main dashboard
            // document.getElementById('viral-sequences').textContent = data.Viral_Hits.length;
            // document.getElementById('known-viruses-count').textContent = knownViruses.length;
        }}

        // Helper function to format E-values nicely
        function formatEvalue(evalue) {{
            if (evalue === null || evalue === undefined) return 'N/A';
            
            // If it's already in scientific notation as a string, return it
            if (typeof evalue === 'string' && evalue.includes('e')) return evalue;
            
            // Convert to number if it's a string
            const evalNum = typeof evalue === 'string' ? parseFloat(evalue) : evalue;
            
            // Format small numbers in scientific notation
            if (evalNum < 0.001) {{
                return evalNum.toExponential(1);
            }}
            
            return evalNum.toString();
        }}

        // Function to scroll to a specific contig when clicking on a virus card
        function scrollToContig(contigId) {{
            document.querySelector(`.visualization-wrapper[data-contig="${{contigId}}"]`)
                .scrollIntoView({{ behavior: 'smooth' }});
        }}


        function toggleIndex() {{
            const index = document.getElementById('contig-index');
            const container = document.querySelector('.page-container');
            index.classList.toggle('collapsed');
            container.classList.toggle('index-collapsed');
        }}


        function createIndex(data) {{
            const indexList = document.getElementById('index-list');
            indexList.innerHTML = ''; // Clear existing items
            
            // Group contigs by sample name
            const sampleGroups = {{}};
            
            // Group viral hits by sample name
            data.Viral_Hits.forEach(hit => {{
                const sampleName = hit.Sample_name || 'Unknown Sample';
                if (!sampleGroups[sampleName]) {{
                    sampleGroups[sampleName] = [];
                }}
                sampleGroups[sampleName].push(hit);
            }});
            
            // Create expandable sections for each sample
            Object.keys(sampleGroups).sort().forEach(sampleName => {{
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
                sampleGroups[sampleName].forEach(hit => {{
                    const li = document.createElement('li');
                    li.className = 'contig-item';
                    
                    // Create a wrapper for the QueryID to allow styling
                    const queryIdSpan = document.createElement('span');
                    queryIdSpan.className = 'query-id';
                    queryIdSpan.textContent = hit.QueryID;
                    
                    // Add indicators for known viruses (high confidence hits)
                    if (hit.vq_score >= 90) {{
                        const indicator = document.createElement('span');
                        indicator.className = 'known-virus-indicator';
                        indicator.title = 'Known virus (Score ≥ 80)';
                        indicator.innerHTML = '★'; // Star icon for known viruses
                        li.appendChild(indicator);
                    }}

                    if (hit.vq_score <= 20) {{
                        const indicator = document.createElement('span');
                        indicator.className = 'no-virus-indicator';
                        // indicator.title = 'Known virus (≥90% identity, ≥70% coverage)';
                        indicator.innerHTML = '✘'; // X icon for no viruses
                        li.appendChild(indicator);
                    }}
                    
                    li.appendChild(queryIdSpan);
                    
                    // Add click event to scroll to the contig visualization
                    li.onclick = () => {{
                        document.querySelector(`.visualization-wrapper[data-contig="${{hit.QueryID}}"]`)
                            .scrollIntoView({{ behavior: 'smooth' }});
                    }};
                    
                    contigList.appendChild(li);
                }});
                
                // Toggle sample group expansion when header is clicked
                sampleHeader.onclick = (e) => {{
                    // Don't trigger if clicking on a child element that has its own click handler
                    if (e.target !== sampleHeader && e.target !== toggleIcon && e.target !== sampleNameElement) {{
                        return;
                    }}
                    
                    sampleGroup.classList.toggle('collapsed');
                    toggleIcon.innerHTML = sampleGroup.classList.contains('collapsed') ? '►' : '▼';
                }};
                
                // Assemble sample group
                sampleGroup.appendChild(sampleHeader);
                sampleGroup.appendChild(contigList);
                indexList.appendChild(sampleGroup);
            }});
        }}

        // Add the necessary CSS for the new index structure
        const style = document.createElement('style');
        style.textContent = `
            .sample-group {{
                margin-bottom: 10px;
            }}
            
            .sample-header {{
                display: flex;
                align-items: center;
                padding: 8px 12px;
                background-color: rgba(0, 0, 0, 0.05);
                border-radius: 4px;
                cursor: pointer;
                font-weight: 600;
                color: var(--text-primary);
                transition: background-color 0.2s;
            }}
            
            .sample-header:hover {{
                background-color: rgba(0, 0, 0, 0.1);
            }}
            
            .toggle-icon {{
                margin-right: 8px;
                font-size: 10px;
                transition: transform 0.2s;
            }}
            
            .sample-name {{
                flex: 1;
            }}
            
            .count-badge {{
                background-color: var(--primary-color);
                color: white;
                border-radius: 12px;
                padding: 2px 8px;
                font-size: 12px;
                font-weight: 500;
            }}
            
            .contig-list {{
                list-style: none;
                padding-left: 20px;
                margin: 5px 0;
                max-height: 400px;
                overflow-y: auto;
            }}
            
            .sample-group.collapsed .contig-list {{
                display: none;
            }}
            
            .contig-item {{
                padding: 6px 10px;
                border-radius: 3px;
                cursor: pointer;
                transition: background-color 0.2s;
                display: flex;
                align-items: center;
            }}
            
            .contig-item:hover {{
                background-color: rgba(0, 0, 0, 0.05);
            }}
            
            .query-id {{
                flex: 1;
            }}
            
            .known-virus-indicator {{
                color: #f39c12;
                margin-right: 6px;
                font-size: 14px;
            }}

            .no-virus-indicator {{
                color: #f32c12;
                margin-right: 6px;
                font-size: 14px;
            }}
        `;
        document.head.appendChild(style);


            function createVisualizations(data) {{
                // container for all visualizations
                const container = d3.select("#visualization-container");
                container.selectAll("*").remove();

                // styling constants
                const styles = {{
                    trackHeight: 130, // accommodate Pfam domains
                    orfHeight: 20,
                    pfamHeight: 20,
                    pfamRadius: 10,
                    margin: {{ top: 50, right: 50, bottom: 100, left: 50 }},
                    colors: {{
                        complete: "#3498db",
                        "5-prime-partial": "#e74c3c",
                        "3-prime-partial": "#f1c40f",
                        genomeLine: "#2c3e50",
                        axis: "#7f8c8d"
                    }},
                    fonts: {{
                        primary: "system-ui, -apple-system, sans-serif",
                        size: {{
                            title: "19px",
                            subtitle: "15px",
                            label: "15px",
                            tooltip: "15px"
                        }}
                    }},
                    pfamColors: [
                        "#ff6b6b", "#58b368", "#bc6ff1", "#ffa048",
                        "#5191d1", "#ffcc29", "#f06595", "#38b2ac"
                    ]
                }};

                    // process each viral hit
                    data.Viral_Hits.forEach((viralHit, index) => {{
                    const contigId = viralHit.QueryID;
                    const contigLength = Math.trunc(viralHit.BLASTx_Qlength || viralHit.BLASTn_Qlength);
                    const orfs = data.ORF_Data[contigId]
                    const fullSeq = viralHit.FullSeq
                    const organism = viralHit.BLASTx_Organism_Name
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
                        .attr("transform", `translate(${{styles.margin.left}}, ${{styles.margin.top}})`);

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
                        .text(`${{contigId}}: ${{contigLength}} nt - ${{organism}}`);

                    // add subtitle
                    const subtitleText = titleGroup.append("text")
                        .attr("class", "visualization-title")
                        .attr("x", 0)
                        .attr("y", 0)
                        .attr("text-anchor", "start")
                        .style("font-family", styles.fonts.primary)
                        .style("font-size", styles.fonts.size.subtitle)
                        .style("font-weight", "300")
                        .text(`Score: ${{vqScore}}`);

                    // bounding box of the title to position the button
                    const titleBBox = titleText.node().getBBox();
                    const subtitleBBox = subtitleText.node().getBBox();
                    const buttonXPosition = titleBBox.width + 10;

                    // copy button
                    const copyButton = titleGroup.append("g")
                        .attr("class", "copy-fasta-button")
                        .style("cursor", "pointer")
                        .on("click", function() {{
                            const fastaHeader = `>${{contigId}}`;
                            const fastaSequence = fullSeq;
                            const fastaString = `${{fastaHeader}}\n${{fastaSequence}}`;

                            navigator.clipboard.writeText(fastaString)
                                .then(() => {{
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
                                }})
                                .catch(err => {{
                                    console.error('Falha ao copiar: ', err);
                                    d3.select(this).select("rect")
                                        .style("fill", "#d9534f");
                                    d3.select(this).select("text")
                                        .style("fill", "white");
                                }});
                        }});

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
                    copyButton.attr("transform", `translate(${{innerWidth - buttonWidth - styles.margin.right + 20}}, -35)`); // margin adjust

                    // adding D3 styles
                    copyButton.on("mouseover", function() {{
                        d3.select(this).select("rect")
                            .transition()
                            .duration(200)
                            .style("fill", hoverBackgroundColor)
                            .style("stroke", hoverBorderColor);
                        d3.select(this).select("text")
                            .transition()
                            .duration(200)
                            .style("fill", hoverTextColor);
                    }})
                    .on("mouseout", function() {{
                        d3.select(this).select("rect")
                            .transition()
                            .duration(200)
                            .style("fill", backgroundColor)
                            .style("stroke", borderColor);
                        d3.select(this).select("text")
                            .transition()
                            .duration(200)
                            .style("fill", textColor);
                    }});


                    // scales
                    const xScale = d3.scaleLinear()
                        .domain([0, parseInt(contigLength)])
                        .range([0, innerWidth]);

                    // track group
                    const trackGroup = mainGroup.append("g");

                    // genome line with gradient
                    const gradientId = `genomeLineGradient-${{contigId}}`;
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
                        .tickFormat(d => `${{d.toLocaleString()}}bp`);

                    trackGroup.append("g")
                        .attr("class", "scale-axis")
                        .attr("transform", `translate(0, ${{styles.trackHeight + 10}})`)
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
                            ${{viralHit.BLASTx_Subject_Title ? `
                                <strong>Subject:</strong> ${{viralHit.BLASTx_Subject_Title}}<br>
                                <strong>Organism:</strong> ${{viralHit.BLASTx_Organism_Name}}<br>
                                <strong>Coverage:</strong> ${{viralHit.BLASTx_Cover}}%<br>
                                <strong>Identity:</strong> ${{viralHit.BLASTx_Ident}}%<br>
                                <strong>E-value:</strong> ${{viralHit.BLASTx_evalue}}<br>
                                <strong>Subj. length:</strong> ${{Math.trunc(viralHit.BLASTx_Slength)}}
                            ` : 'No BLASTx hits'}}
                        `);

                    // BLASTn column
                    blastInfo.append("div")
                        .html(`
                            <strong style="color: #2c3e50; font-size: 14px">BLASTn Results</strong>
                            <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                            ${{viralHit.BLASTn_Subject_Title ? `
                                <strong>Subject:</strong> ${{viralHit.BLASTn_Subject_Title}}<br>
                                <strong>Coverage:</strong> ${{viralHit.BLASTn_Cover}}%<br>
                                <strong>Identity:</strong> ${{viralHit.BLASTn_Ident}}%<br>
                                <strong>E-value:</strong> ${{viralHit.BLASTn_evalue}}<br>
                                <strong>Subj. length:</strong> ${{Math.trunc(viralHit.BLASTn_Slength)}}
                            ` : 'No BLASTn hits'}}
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
                                <strong style="color: #555; font-weight: 600; margin-right: 5px;">Scientific Name:</strong> ${{viralHit.ScientificName}}
                            </div>
                            <div style="margin-bottom: 8px;">
                                <strong style="color: #555; font-weight: 600; margin-right: 5px;">Species:</strong> <span style="color: ${{viralHit.Species ? '#333' : '#999'}}; font-style: ${{viralHit.Species ? 'normal' : 'italic'}};">${{viralHit.Species || 'N/A'}}</span>
                            </div>
                            <div style="margin-bottom: 8px;">
                                <strong style="color: #555; font-weight: 600; margin-right: 5px;">Rank:</strong> <span style="color: ${{viralHit.NoRank ? '#333' : '#999'}}; font-style: ${{viralHit.NoRank ? 'normal' : 'italic'}};">${{viralHit.NoRank || 'N/A'}}</span>
                            </div>
                            <br>
                            <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 10px;">                           
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">Kingdom</strong>
                                    <span style="color: ${{viralHit.Kingdom ? '#333' : '#999'}}; font-style: ${{viralHit.Kingdom ? 'normal' : 'italic'}}; display: block;">${{viralHit.Kingdom || 'N/A'}}</span>
                                </div>
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">Phylum</strong>
                                    <span style="color: ${{viralHit.Phylum ? '#333' : '#999'}}; font-style: ${{viralHit.Phylum ? 'normal' : 'italic'}}; display: block;">${{viralHit.Phylum || 'N/A'}}</span>
                                </div>
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">Class</strong>
                                    <span style="color: ${{viralHit.Class ? '#333' : '#999'}}; font-style: ${{viralHit.Class ? 'normal' : 'italic'}}; display: block;">${{viralHit.Class || 'N/A'}}</span>
                                </div>
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">Order</strong>
                                    <span style="color: ${{viralHit.Order ? '#333' : '#999'}}; font-style: ${{viralHit.Order ? 'normal' : 'italic'}}; display: block;">${{viralHit.Order || 'N/A'}}</span>
                                </div>
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">Family</strong>
                                    <span style="color: ${{viralHit.Family ? '#333' : '#999'}}; font-style: ${{viralHit.Family ? 'normal' : 'italic'}}; display: block;">${{viralHit.Family || 'N/A'}}</span>
                                </div>
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">Genome</strong>
                                    <span style="color: ${{viralHit.Genome ? '#333' : '#999'}}; font-style: ${{viralHit.Genome ? 'normal' : 'italic'}}; display: block;">${{viralHit.Genome || 'N/A'}}</span>
                                </div>
                                <div>
                                    <strong style="color: #555; font-weight: 600; display: block;">TaxID</strong>
                                    <span style="color: ${{viralHit.TaxId ? '#333' : '#999'}}; font-style: ${{viralHit.TaxId ? 'normal' : 'italic'}}; display: block;">${{viralHit.TaxId ? Math.trunc(viralHit.TaxId) : 'N/A'}}</span>
                                </div>                         
                            </div>
                        `);

                    // Pfam domain summary panel
                    const hmmHits = data.HMM_hits.filter(hit => hit.QueryID === contigId);
                    if (hmmHits.length > 0) {{
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
                        hmmHits.forEach(hit => {{
                            // extract keys that match Pfam_Accession
                            const pfamKeys = Object.keys(hit).filter(key => key.match(/^Pfam_Accession(_\d+)?$/));

                            // Pfam domain information
                            pfamKeys.forEach(accKey => {{
                                const suffix = accKey.replace('Pfam_Accession', '');
                                if (hit[accKey] && hit[`Pfam_Description${{suffix}}`]) {{
                                    pfamDomains.push({{
                                        accession: hit[accKey],
                                        description: hit[`Pfam_Description${{suffix}}`],
                                        type: hit[`Pfam_Type${{suffix}}`] || 'Unknown',
                                        orf: hit.Query_name
                                    }});
                                }}
                            }});
                        }});

                        // remove duplicates and show unique Pfam domain
                        const uniquePfamDomains = Array.from(new Map(
                            pfamDomains.map(domain => [domain.accession, domain])
                        ).values());

                        uniquePfamDomains.forEach((domain, i) => {{
                            const colorIndex = i % styles.pfamColors.length;
                            const domainColor = styles.pfamColors[colorIndex];

                            pfamList.append("div")
                                .style("padding", "5px")
                                .style("border-left", `4px solid ${{domainColor}}`)
                                .style("background-color", "#fff")
                                .html(`
                                    <strong>${{domain.accession}}</strong>: ${{domain.description}}<br>
                                    <small>Type: ${{domain.type}}</small>
                                `);
                        }});
                    }}

                    if (viralHit.AI_summary) {{
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
                    }}



                    // tooltip (one per visualization)
                    const tooltip = d3.select("body").append("div")
                        .attr("class", `tooltip-${{contigId}}`)
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
                    const createTooltipContent = (orf, hmmHit) => {{
                        const sequence = orf.sequence || '';
                        const truncatedSeq = sequence.substring(0, 8) + "..";
                        const buttonId = `copy-button-${{orf.Query_name.replace(/[^a-zA-Z0-9]/g, '-')}}`;
                        const feedbackId = `feedback-${{orf.Query_name.replace(/[^a-zA-Z0-9]/g, '-')}}`;


                        // position string based on strand
                        const positionString = orf.strand === '+' ? 
                            `${{orf.start}}-${{orf.end}}` : 
                            `${{orf.end}}-${{orf.start}}`;

                        return `
                            <div style="font-family: ${{styles.fonts.primary}}">
                                <strong style="color: #2c3e50">${{orf.Query_name}}</strong><br>
                                <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                                <strong>Position (nt):</strong> ${{positionString}}<br>
                                <strong>Type:</strong> ${{orf.type}}<br>
                                <strong>Length (nt):</strong> ${{orf.length}}<br>
                                <strong>Length (aa):</strong> ${{orf.length / 3}}<br>
                                <strong>Frame:</strong> ${{orf.frame}}<br>
                                <strong>Strand:</strong> ${{orf.strand}}<br>
                                <strong>Codons:</strong> ${{orf.start_codon}} → ${{orf.stop_codon}}<br>
                                <strong>Sequence (aa):</strong>
                                <span style="font-family: monospace">${{truncatedSeq}}</span>
                                <button
                                    id="${{buttonId}}"
                                    class="copy-button"
                                    data-sequence="${{sequence}}"
                                    onclick="handleCopyClick(this)">
                                    Copy
                                </button>
                                <span id="${{feedbackId}}" class="copy-feedback">
                                    Copied!
                                </span>
                                ${{hmmHit ? `
                                    <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                                    <strong style="color: #2c3e50">HMM Hits:</strong><br>
                                    ${{['RVDB', 'Vfam', 'EggNOG'].map(db =>
                                        hmmHit[`${{db}}_TargetID`] ?
                                            `<strong>${{db}}:</strong> ${{hmmHit[`${{db}}_TargetID`]}}
                                            ${{hmmHit[`${{db}}_Description`] ? ` - ${{hmmHit[`${{db}}_Description`]}}` : ''}}<br>`
                                        : ''
                                    ).join('')}}
                                ` : ''}}
                            </div>
                        `;
                    }};

                    // tooltip content for Pfam domains
                    const createPfamTooltipContent = (pfamData) => {{
                        return `
                            <div class="pfam-tooltip" style="font-family: ${{styles.fonts.primary}}">
                                <strong style="color: #2c3e50">${{pfamData.Pfam_Accession}}</strong><br>
                                <hr style="border: 1px solid #ecf0f1; margin: 5px 0">
                                <strong>Description:</strong> ${{pfamData.Pfam_Description}}<br>
                                <strong>Type:</strong> ${{pfamData.Pfam_Type}}<br>
                                <strong>Position (aa):</strong> ${{pfamData.Pfam_Start}}-${{pfamData.Pfam_End}}<br>
                                <strong>Length (aa):</strong> ${{pfamData.Pfam_length}}<br>
                                <strong>Score:</strong> ${{parseFloat(pfamData.Pfam_Score).toFixed(2)}}<br>
                                ${{pfamData.Pfam_Info ? `
                                    <div class="pfam-info">
                                        ${{pfamData.Pfam_Info.replace(/<\/?p>/g, '')}}
                                    </div>
                                ` : ''}}
                            </div>
                        `;
                    }};

                    // handle copy button click
                    window.handleCopyClick = function(button) {{
                        const sequence = button.getAttribute('data-sequence');
                        const feedbackId = button.id.replace('copy-button-', 'feedback-');
                        const feedbackElement = document.getElementById(feedbackId);

                        navigator.clipboard.writeText(sequence)
                            .then(() => {{
                                // feedback
                                button.style.display = 'none';
                                feedbackElement.style.display = 'inline';
                                feedbackElement.classList.add('show');

                                // reset after 2 seconds
                                setTimeout(() => {{
                                    button.style.display = 'inline';
                                    feedbackElement.classList.remove('show');
                                    setTimeout(() => {{
                                        feedbackElement.style.display = 'none';
                                    }}, 300);
                                }}, 2000);
                            }})
                            .catch(err => {{
                                console.error('Failed to copy:', err);
                                feedbackElement.textContent = 'Error copying!';
                                feedbackElement.style.color = '#e74c3c';
                                feedbackElement.style.display = 'inline';
                                feedbackElement.classList.add('show');
                            }});
                    }};

                    // show tooltip
                    const showTooltip = (event, content) => {{
                        tooltip
                            .style("display", "block")
                            .style("left", (event.pageX + 15) + "px")
                            .style("top", (event.pageY - 10) + "px")
                            .html(content);
                    }};

                    // hide tooltip
                    const hideTooltip = () => {{
                        if (!isTooltipPinned) {{
                            tooltip.style("display", "none");
                            activeTooltip = null;
                        }}
                    }};

                    // process ORFs
                    orfs.forEach(orf => {{
                        const orfGroup = trackGroup.append("g");
                        const yOffset = orf.frame > 0 ? -0 : 0;
                        const orfWidth = xScale(orf.end) - xScale(orf.start);
                        const arrowSize = Math.min(15, orfWidth / 5);

                        // arrow-shaped ORF path
                        const orfPath = orf.strand === "+" ?
                            `M ${{xScale(orf.start)}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset}}
                            L ${{xScale(orf.end) - arrowSize}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset}}
                            L ${{xScale(orf.end)}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight / 2}}
                            L ${{xScale(orf.end) - arrowSize}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}}
                            L ${{xScale(orf.start)}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}}
                            Z` :
                            `M ${{xScale(orf.start)}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight / 2}}
                            L ${{xScale(orf.start) + arrowSize}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset}}
                            L ${{xScale(orf.end)}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset}}
                            L ${{xScale(orf.end)}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}}
                            L ${{xScale(orf.start) + arrowSize}} ${{(styles.trackHeight - styles.orfHeight) / 2 + yOffset + styles.orfHeight}}
                            Z`;

                        // ORF shape with interaction
                        orfGroup.append("path")
                            .attr("d", orfPath)
                            .style("fill", styles.colors[orf.type])
                            .style("filter", "none")
                            .style("stroke", "#333")
                            .style("cursor", "pointer")
                            .on("mouseover", (event) => {{
                                const hmmHit = data.HMM_hits.find(hit => hit.Query_name === orf.Query_name);
                                isTooltipPinned = false;
                                showTooltip(event, createTooltipContent(orf, hmmHit));
                            }})
                            .on("mousemove", (event) => {{
                                if (!isTooltipPinned) {{
                                    tooltip
                                        .style("left", (event.pageX + 15) + "px")
                                        .style("top", (event.pageY - 10) + "px");
                                }}
                            }})
                            .on("mouseout", () => {{
                                tooltipTimeout = setTimeout(hideTooltip, 3000);
                            }})
                            .on("click", (event) => {{
                                const hmmHit = data.HMM_hits.find(hit => hit.Query_name === orf.Query_name);
                                isTooltipPinned = !isTooltipPinned;
                                if (isTooltipPinned) {{
                                    showTooltip(event, createTooltipContent(orf, hmmHit));
                                    if (tooltipTimeout) {{
                                        clearTimeout(tooltipTimeout);
                                    }}
                                }} else {{
                                    hideTooltip();
                                }}
                                event.stopPropagation();
                            }});

                        // corresponding HMM hit and extract Pfam domains
                        const hmmHit = data.HMM_hits.find(hit => hit.Query_name === orf.Query_name);
                        if (hmmHit) {{
                            // all Pfam domains from the HMM hit
                            const pfamDomains = [];

                            // keys that match Pfam_Accession pattern
                            const pfamKeys = Object.keys(hmmHit).filter(key => key.match(/^Pfam_Accession(_\d+)?$/));

                            // domain object with all related information
                            pfamKeys.forEach(accKey => {{
                                // suffix (empty string or _2, _3, etc.)
                                const suffix = accKey.replace('Pfam_Accession', '');

                                // process if there's a start and end position
                                if (hmmHit[`Pfam_Start${{suffix}}`] && hmmHit[`Pfam_End${{suffix}}`]) {{
                                    pfamDomains.push({{
                                        Pfam_Accession: hmmHit[accKey],
                                        Pfam_Description: hmmHit[`Pfam_Description${{suffix}}`] || '',
                                        Pfam_Info: hmmHit[`Pfam_Info${{suffix}}`] || '',
                                        Pfam_Type: hmmHit[`Pfam_Type${{suffix}}`] || '',
                                        Pfam_Score: hmmHit[`Pfam_Score${{suffix}}`] || '0',
                                        Pfam_Start: parseFloat(hmmHit[`Pfam_Start${{suffix}}`]),
                                        Pfam_End: parseFloat(hmmHit[`Pfam_End${{suffix}}`]),
                                        Pfam_length: parseFloat(hmmHit[`Pfam_length${{suffix}}`] || '0')
                                    }});
                                }}
                            }});

                            // Pfam domains below the ORF
                            const pfamGroup = trackGroup.append("g");

                            pfamDomains.forEach((domain, i) => {{
                                const colorIndex = i % styles.pfamColors.length;
                                const domainColor = styles.pfamColors[colorIndex];

                                // Pfam domain positions from amino acid to nucleotide coordinates
                                let domainStartNuc, domainEndNuc;

                                // frame offset (0, 1, or 2 depending on the reading frame)
                                const frameOffset = Math.abs(orf.frame) - 1; // Convert frame (1-6) to offset (0-2)

                                if (orf.strand === "+") {{
                                    // for positive strand, translate directly from AA to nucl
                                    // each amino acid is 3 nucleotides
                                    domainStartNuc = orf.start + frameOffset + (domain.Pfam_Start - 1) * 3;
                                    domainEndNuc = orf.start + frameOffset + (domain.Pfam_End * 3) - 1; // -1 because Pfam_End is inclusive
                                }} else {{
                                    // for strand == -, we count from the end backwards
                                    // We need to reverse the amino acid coordinates
                                    domainStartNuc = orf.end - frameOffset - (domain.Pfam_End * 3) + 1; // +1 because we're counting backwards
                                    domainEndNuc = orf.end - frameOffset - (domain.Pfam_Start - 1) * 3;
                                }}

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
                                    .on("mouseover", (event) => {{
                                        showTooltip(event, createPfamTooltipContent(domain));
                                    }})
                                    .on("mousemove", (event) => {{
                                        if (!isTooltipPinned) {{
                                            tooltip
                                                .style("left", (event.pageX + 15) + "px")
                                                .style("top", (event.pageY - 10) + "px");
                                        }}
                                    }})
                                    .on("mouseout", () => {{
                                        tooltipTimeout = setTimeout(hideTooltip, 3000);
                                    }})
                                    .on("click", (event) => {{
                                        isTooltipPinned = !isTooltipPinned;
                                        if (isTooltipPinned) {{
                                            showTooltip(event, createPfamTooltipContent(domain));
                                            if (tooltipTimeout) {{
                                                clearTimeout(tooltipTimeout);
                                            }}
                                        }} else {{
                                            hideTooltip();
                                        }}
                                        event.stopPropagation();
                                    }});

                                // domain label if space allows
                                if (xScale(domainEndNuc) - xScale(domainStartNuc) > 80) {{
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
                                }}
                            }});
                        }}
                    }});
                }});

                // hide loading indicator
                document.getElementById("loading").style.display = "none";

                // add event listener to close tooltip when clicking elsewhere
                document.addEventListener("click", (event) => {{
                    if (!event.target.closest(".tooltip")) {{
                        d3.selectAll("[class^='tooltip-']").style("display", "none");
                        isTooltipPinned = false;
                    }}
                }});
            }}

            // load data and initialize visualization
            document.addEventListener("DOMContentLoaded", () => {{

            // populate the dashboard
            function populateDashboard(data) {{
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
                document.getElementById('blastn-db').textContent = `BLASTn: ${{data.args.blastn}}`;
                document.getElementById('blastx-db').textContent = `BLASTx: ${{data.args.diamond_blastx}}`;
                document.getElementById('rvdb_hmm').textContent = `RVDB: ${{data.args.rvdb_hmm}}`;
                document.getElementById('pfam_hmm').textContent = `Pfam: ${{data.args.pfam_hmm}}`;

                // pie chart for sequence distribution
                createSequenceDistributionChart(data);
            }}

            // create sequence distribution chart
            function createSequenceDistributionChart(data) {{
                const width = document.getElementById('sequences-chart').clientWidth;
                const height = 200;
                const radius = Math.min(width, height) / 2 * 0.8;

                const svg = d3.select("#sequences-chart")
                    .append("svg")
                    .attr("width", width)
                    .attr("height", height)
                    .append("g")
                    .attr("transform", `translate(${{width/2}}, ${{height/2}})`);

                // data for pie chart
                const pieData = [
                    {{ label: "Viral", value: data.number_viralSeqs, color: "#e74c3c" }},
                    {{ label: "Non-viral ≥ 500nt", value: data.number_filteredSeqs - data.number_viralSeqs, color: "#3498db" }},
                    {{ label: "Contigs < 500nt", value: data.number_originalSeqs - data.number_filteredSeqs, color: "#95a5a6" }}
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
                    .on("mouseover", function() {{
                        d3.select(this)
                            .style("opacity", 1);
                    }})
                    .on("mouseout", function() {{
                        d3.select(this)
                            .style("opacity", 0.8);
                    }});

                    // background rectangles for percentage labels
                    svg.selectAll("rect.label-bg") // Seleciona elementos rect com a classe label-bg (se existirem)
                        .data(pie(pieData))
                        .enter()
                        .append("rect")
                        .attr("class", "label-bg") // add class to select style
                        .attr("transform", d => `translate(${{arcLabels.centroid(d)}})`) // arc centroid
                        .attr("width", d => {{ // width based on text
                            const percent = (d.data.value / data.number_originalSeqs * 100).toFixed(1);
                            const textLength = percent > 0.00000001 ? `${{percent}}%`.length * 8 : 0; // width estimative 
                            return textLength + 10; // add a few padding
                        }})
                        .attr("height", 16) // fix height 
                        .attr("x", d => {{ // adjust x to centralize rectangle
                            const percent = (d.data.value / data.number_originalSeqs * 100).toFixed(1);
                            const textLength = percent > 0.00000001 ? `${{percent}}%`.length * 8 : 0;
                            return - (textLength + 10) / 2;
                        }})
                        .attr("y", -8) // adjust y to centralize rectangle
                        .attr("fill", "white")
                        .attr("rx", 5) // radius for round corners
                        .attr("ry", 5);

                // percentage labels
                svg.selectAll("text")
                    .data(pie(pieData))
                    .enter()
                    .append("text")
                    .attr("transform", d => `translate(${{arcLabels.centroid(d)}})`)
                    .attr("text-anchor", "middle")
                    .attr("dominant-baseline", "middle")
                    .text(d => {{
                        const percent = (d.data.value / data.number_originalSeqs * 100).toFixed(1);
                        return percent > 0.00000001 ? `${{percent}}%` : '';
                    }})
                    .style("font-size", "12px")
                    .style("fill", "black")
                    .style("font-weight", "bold");

                // legend
                const legend = svg.selectAll(".legend")
                    .data(pieData)
                    .enter()
                    .append("g")
                    .attr("class", "legend")
                    .attr("transform", (d, i) => `translate(${{width/2 - 180}}, ${{i * 20 - 30}})`);

                legend.append("rect")
                    .attr("width", 12)
                    .attr("height", 12)
                    .attr("fill", d => d.color);

                legend.append("text")
                    .attr("x", 20)
                    .attr("y", 6)
                    .attr("dy", ".35em")
                    .style("font-size", "12px")
                    .text(d => `${{d.label}} (${{d.value}})`);
            }}





            // graph plot
            function createTaxonomyVisualization(data) {{
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
                    .attr("viewBox", `0 0 ${{width}} ${{height-210}}`)
                    .attr("preserveAspectRatio", "xMidYMid meet");
                
                // Main group
                const mainGroup = svg.append("g")
                    .attr("transform", `translate(${{width / 2}}, ${{height / 3.4}})`);
                
                // Enhanced styling constants for more compact visualization
                const styles = {{
                    margin: {{ top: 100, right: 10, bottom: 10, left: 10 }},
                    colors: {{
                        phylum: "#3498db",     // Blue for phylum
                        family: "#2ecc71",     // Green for family
                        virus: "#e74c3c",      // Red for viruses
                        link: "#ecf0f1",       // Light gray for links
                        highlight: "#ff7f0e"   // Orange for highlighting
                    }},
                    fonts: {{
                        primary: "system-ui, -apple-system, sans-serif",
                        size: {{
                            title: "15px",     // Slightly smaller
                            label: "11px",     // Slightly smaller
                            tooltip: "13px"    // Slightly smaller
                        }}
                    }},
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
                    force: {{
                        center: 0.1,         // Strength of center force
                        charge: 90,         // Base charge strength
                        link: 30,            // Base link distance
                        boundary: 0.15,      // Strength of boundary force
                        radius: 30         // Boundary radius
                    }}
                }};

                // Process the data to extract unique phyla, families, and viruses
                // This is a placeholder; you'll need to adapt this to your actual data structure
                const processedData = processViralData(data);
                
                // Create force simulation with more constraints for compact layout
                const simulation = d3.forceSimulation(processedData.nodes)
                    .force("link", d3.forceLink(processedData.links).id(d => d.id).distance(d => {{
                        // Shorter distances for a more compact layout
                        if (d.source.type === "phylum" && d.target.type === "family") return 220;
                        return 90;
                    }}))
                    .force("charge", d3.forceManyBody().strength(d => {{
                        // Less repulsion for more compact layout
                        if (d.type === "phylum") return -250;
                        if (d.type === "family") return -250;
                        return -70;
                    }}))
                    .force("center", d3.forceCenter(0, 0))
                    // Add force to keep nodes within bounds
                    .force("x", d3.forceX().strength(0.1))
                    .force("y", d3.forceY().strength(0.1))
                    .force("collision", d3.forceCollide().radius(d => {{
                        // Smaller collision radiuses for compact layout
                        if (d.type === "phylum") return 50;
                        if (d.type === "family") return 50;
                        return 20;
                    }}));

                // Create links
                const link = mainGroup.append("g")
                    .attr("class", "links")
                    .selectAll("line")
                    .data(processedData.links)
                    .enter().append("line")
                    .attr("stroke", styles.colors.link)
                    .attr("stroke-opacity", 0.6)
                    .attr("stroke-width", d => {{
                        // Different stroke widths based on link type
                        if (d.source.type === "phylum" && d.target.type === "family") return 2;
                        return 1;
                    }});

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
                    .attr("r", d => {{
                        // Smaller sizes for compact layout
                        if (d.type === "phylum") return 20;
                        if (d.type === "family") return 12;
                        return 5;
                    }})
                    .attr("fill", d => {{
                        // Different colors based on node type
                        if (d.type === "phylum") return styles.phylumColorScale(d.id);
                        if (d.type === "family") return styles.familyColorScale(d.id);
                        return styles.colors.virus;
                    }})
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 1.5);

                // Add labels only to phylum and family nodes (omit viral names)
                node.filter(d => d.type !== "virus") // Only add text to phylum and family nodes
                    .append("text")
                    .text(d => d.name)
                    .attr("x", 0)
                    .attr("y", d => {{
                        // Position text differently based on node type
                        if (d.type === "phylum") return -25;
                        if (d.type === "family") return -15;
                        return 0;
                    }})
                    .attr("text-anchor", "middle")
                    .attr("font-family", styles.fonts.primary)
                    .attr("font-size", d => {{
                        // Different font sizes based on node type
                        if (d.type === "phylum") return styles.fonts.size.title;
                        return styles.fonts.size.label;
                    }})
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
                node.on("mouseover", function(event, d) {{
                    d3.select(this).select("circle")
                        .transition()
                        .duration(200)
                        .attr("r", d => {{
                            if (d.type === "phylum") return 25;
                            if (d.type === "family") return 16;
                            return 8; // Slightly larger for viruses but still compact
                        }});

                    // Show tooltip with information - show viral names here instead of on the graph
                    tooltip.style("visibility", "visible")
                        .html(() => {{
                            let content = `<strong>${{d.name}}</strong><br>Type: ${{d.type.charAt(0).toUpperCase() + d.type.slice(1)}}`;
                            if (d.type === "virus") {{
                                content += `<br>Family: ${{d.family || 'N/A'}}`;
                                content += `<br>Phylum: ${{d.phylum || 'N/A'}}`;
                                // Add more virus details if available
                                if (d.Genome) content += `<br>Genome Type: ${{d.Genome}}`;
                                if (d.taxId) content += `<br>TaxID: ${{d.taxId}}`;
                            }} else if (d.type === "family") {{
                                content += `<br>Phylum: ${{d.phylum || 'N/A'}}`;
                                // Add count of viruses in this family
                                const virusCount = processedData.links.filter(link => 
                                    link.source.id === d.id || 
                                    (typeof link.source === 'string' && link.source === d.id)
                                ).length;
                                content += `<br>Viruses: ${{virusCount}}`;
                            }} else if (d.type === "phylum") {{
                                // Add count of families in this phylum
                                const familyCount = processedData.links.filter(link => 
                                    (link.source.type === "phylum" && link.source.id === d.id) || 
                                    (typeof link.source === 'string' && link.source === d.id && link.target.type === "family")
                                ).length;
                                content += `<br>Families: ${{familyCount}}`;
                            }}
                            return content;
                        }});

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
                }})
                .on("mousemove", function(event) {{
                    tooltip.style("top", (event.pageY - 10) + "px")
                        .style("left", (event.pageX + 10) + "px");
                }})
                .on("mouseout", function(event, d) {{
                    d3.select(this).select("circle")
                        .transition()
                        .duration(200)
                        .attr("r", d => {{
                            if (d.type === "phylum") return 20;
                            if (d.type === "family") return 12;
                            return 5;
                        }});
                    
                    tooltip.style("visibility", "hidden");

                    // Reset all links and nodes
                    link.style("opacity", 0.6)
                        .style("stroke-width", d => d.source.type === "phylum" ? 2 : 1)
                        .style("stroke", styles.colors.link);
                    
                    node.style("opacity", 1);
                    
                    tooltip.style("visibility", "hidden");

                }});



                // constrain the graph in  a rectangular area
                const bounds = {{
                    x: {{ min: -width/2 + 40, max: width/2 - 40 }},
                    y: {{ min: -height/2 + 40, max: height/2 - 20 }}
                }};

                simulation.on("tick", () => {{
                    // Constrain nodes to a rectangular area to prevent flying away
                    processedData.nodes.forEach(d => {{
                        // Apply rectangular constraints
                        d.x = Math.max(bounds.x.min, Math.min(bounds.x.max, d.x));
                        d.y = Math.max(bounds.y.min, Math.min(bounds.y.max, d.y));
                    }});
                    
                    link
                        .attr("x1", d => d.source.x)
                        .attr("y1", d => d.source.y)
                        .attr("x2", d => d.target.x)
                        .attr("y2", d => d.target.y);

                    node.attr("transform", d => `translate(${{d.x}},${{d.y}})`);
                }});

                // Drag functions
                function dragstarted(event, d) {{
                    if (!event.active) simulation.alphaTarget(0.3).restart();
                    d.fx = d.x;
                    d.fy = d.y;
                }}

                function dragged(event, d) {{
                    d.fx = event.x;
                    d.fy = event.y;
                }}

                function dragended(event, d) {{
                    if (!event.active) simulation.alphaTarget(0);
                    d.fx = null;
                    d.fy = null;
                }}

                // Function to process viral data
                function processViralData(data) {{
                    const nodes = [];
                    const links = [];
                    const phylaMap = new Map();
                    const familyMap = new Map();
                    const virusMap = new Map();  // Track unique viruses to avoid duplicates
                    
                    // Count objects to size nodes proportionally
                    const phylumCounts = {{}};
                    const familyCounts = {{}};
                    
                    // First pass: count items for sizing
                    if (data && data.Viral_Hits) {{
                        data.Viral_Hits.forEach(viralHit => {{
                            const phylum = viralHit.Phylum || "Unknown Phylum";
                            const family = viralHit.Family || "Unknown Family";
                            
                            // Count phyla and families for sizing nodes
                            phylumCounts[phylum] = (phylumCounts[phylum] || 0) + 1;
                            
                            const familyKey = `${{family}}_${{phylum}}`;
                            familyCounts[familyKey] = (familyCounts[familyKey] || 0) + 1;
                        }});
                    }}
                    
                    // Second pass: create nodes and links
                    if (data && data.Viral_Hits) {{
                        data.Viral_Hits.forEach(viralHit => {{
                            const virusName = viralHit.ScientificName || viralHit.QueryID;
                            const phylum = viralHit.Phylum || "Unknown Phylum";
                            const family = viralHit.Family || "Unknown Family";
                            
                            // Add phylum if not exists
                            if (!phylaMap.has(phylum)) {{
                                phylaMap.set(phylum, true);
                                nodes.push({{
                                    id: phylum,
                                    name: phylum,
                                    type: "phylum",
                                    count: phylumCounts[phylum] || 1,
                                    viralHits: [] // To store refs to viral hits in this phylum
                                }});
                            }}
                            
                            // Add family if not exists
                            const familyId = `${{family}}_${{phylum}}`;
                            if (!familyMap.has(familyId)) {{
                                familyMap.set(familyId, true);
                                nodes.push({{
                                    id: familyId,
                                    name: family,
                                    type: "family",
                                    phylum: phylum,
                                    count: familyCounts[familyId] || 1,
                                    viralHits: [] // To store refs to viral hits in this family
                                }});
                                
                                // Link family to phylum
                                links.push({{
                                    source: phylum,
                                    target: familyId,
                                    value: familyCounts[familyId] || 1 // Width based on count
                                }});
                            }}
                            
                            // Add virus (only if unique)
                            const virusId = `${{virusName}}_${{familyId}}`;
                            if (!virusMap.has(virusId)) {{
                                virusMap.set(virusId, true);
                                nodes.push({{
                                    id: virusId,
                                    name: virusName,
                                    type: "virus",
                                    family: family,
                                    phylum: phylum,
                                    Genome: viralHit.Genome,
                                    taxId: viralHit.TaxId ? Math.trunc(viralHit.TaxId) : null,
                                    viralHit: viralHit // Store reference to original data
                                }});
                                
                                // Find the family node and add this virus to its viralHits
                                const familyNode = nodes.find(node => node.id === familyId);
                                if (familyNode) {{
                                    familyNode.viralHits.push(viralHit);
                                }}
                                
                                // Find the phylum node and add this virus to its viralHits
                                const phylumNode = nodes.find(node => node.id === phylum);
                                if (phylumNode) {{
                                    phylumNode.viralHits.push(viralHit);
                                }}
                                
                                // Link virus to family
                                links.push({{
                                    source: familyId,
                                    target: virusId,
                                    value: 1
                                }});
                            }}
                        }});
                    }}
                    
                    return {{ nodes, links }};
                }}
                
                // Add search functionality
                // Add filter controls instead of just search
                const controlPanel = taxonomyDiv.append("div")
                    .style("position", "absolute")
                    .style("top", "30px")
                    .style("left", "71.8%")
                    .style("z-index", "100")
                    .style("background-color", "rgba(255, 255, 255, 0.8)")
                    .style("padding", "10px")
                    .style("border-radius", "4px")
                    .style("box-shadow", "0 1px 3px rgba(0,0,0,0.12)");
                
                // Search box for all node types
                controlPanel.append("input")
                    .attr("type", "text")
                    .attr("placeholder", "Search taxonomy...")
                    .style("padding", "6px 10px")
                    .style("border", "1px solid #ddd")
                    .style("border-radius", "4px")
                    .style("width", "180px")
                    .style("font-family", styles.fonts.primary)
                    .style("font-size", "13px")
                    .style("margin-bottom", "8px")
                    .style("display", "block")
                    .on("input", function() {{
                        const searchTerm = this.value.toLowerCase();
                        
                        // Highlight matching nodes
                        node.each(function(d) {{
                            const element = d3.select(this);
                            if (d.name.toLowerCase().includes(searchTerm)) {{
                                element.select("circle")
                                    .transition()
                                    .duration(200)
                                    .attr("r", d => {{
                                        if (d.type === "phylum") return 25;
                                        if (d.type === "family") return 16;
                                        return 8;
                                    }})
                                    .attr("stroke", "#3498db")
                                    .attr("stroke-width", 2);
                                
                                // Make connections to highlighted nodes more visible
                                link.filter(l => l.source.id === d.id || l.target.id === d.id)
                                    .attr("stroke", "#3498db")
                                    .attr("stroke-opacity", 0.8)
                                    .attr("stroke-width", 1.5);
                            }} else {{
                                element.select("circle")
                                    .transition()
                                    .duration(200)
                                    .attr("r", d => {{
                                        if (d.type === "phylum") return 20;
                                        if (d.type === "family") return 12;
                                        return 5;
                                    }})
                                    .attr("stroke", "#fff")
                                    .attr("stroke-width", 1.5);
                                
                                // Reset link appearance
                                link.filter(l => l.source.id === d.id || l.target.id === d.id)
                                    .attr("stroke", styles.colors.link)
                                    .attr("stroke-opacity", 0.6)
                                    .attr("stroke-width", d => {{
                                        if (d.source.type === "phylum" && d.target.type === "family") return 1.5;
                                        return 1;
                                    }});
                            }}
                        }});
                    }});
                
                // Reset button
                controlPanel.append("button")
                    .text("Reset View")
                    .style("padding", "6px 10px")
                    .style("background-color", "#f5f5f5")
                    .style("border", "1px solid #ddd")
                    .style("border-radius", "4px")
                    .style("cursor", "pointer")
                    .style("font-family", styles.fonts.primary)
                    .style("font-size", "13px")
                    .style("width", "200px")
                    .style("display", "block")
                    .on("click", function() {{
                        // Reset node appearance
                        node.select("circle")
                            .transition()
                            .duration(200)
                            .attr("r", d => {{
                                if (d.type === "phylum") return 20;
                                if (d.type === "family") return 12;
                                return 5;
                            }})
                            .attr("stroke", "#fff")
                            .attr("stroke-width", 1.5);
                        
                        // Reset link appearance
                        link.attr("stroke", styles.colors.link)
                            .attr("stroke-opacity", 0.6)
                            .attr("stroke-width", d => {{
                                if (d.source.type === "phylum" && d.target.type === "family") return 1.5;
                                return 1;
                            }});
                        
                        // Reset simulation
                        simulation.alpha(0.3).restart();
                        
                        // Clear search input
                        controlPanel.select("input").property("value", "");
                    }});
            }}



            // call function
            //populateDashboard(logData);
            populateKnownVirusesDashboard(jsonData);
            createTaxonomyVisualization(jsonData);


                // dealing with JSON dataase
                try {{
                    createIndex(jsonData);
                    createVisualizations(jsonData);
                }} catch (error) {{
                    console.error("Error loading visualizations:", error);
                    document.getElementById("loading").style.display = "none";
                    document.getElementById("error").style.display = "block";
                }}
            }});


                function initializeSelectionAndExport() {{
                    // CSS styles for new elements
                    const styleEl = document.createElement('style');
                    styleEl.textContent = `
                        .checkbox-container {{
                            position: absolute; 
                            left: 14px;
                            top: 14px;
                            z-index: 10;
                        }}
                        .sequence-checkbox {{
                            width: 20px;
                            height: 20px;
                            cursor: pointer;
                        }}
                        .export-panel {{
                            margin: 20px 0;
                            padding: 15px;
                            background-color: #ffffff;
                            border: 1px solid #ecf0f1;
                            border-radius: 10px;
                            display: flex;
                            justify-content: space-between;
                            box-shadow: 0 2px 4px var(--shadow-color);
                        }}
                        .export-btn {{
                            background-color: #3498db;
                            color: white;
                            border: none;
                            padding: 10px 15px;
                            border-radius: 4px;
                            cursor: pointer;
                            font-family: system-ui, -apple-system, sans-serif;
                            font-size: 16px;
                        }}
                        .export-btn:hover {{
                            background-color: #2980b9;
                        }}
                        .select-btns {{
                            display: flex;
                            gap: 10px;
                        }}
                        .select-btn {{
                            background-color: #f1f1f1;
                            border: 1px solid #ddd;
                            padding: 8px 12px;
                            border-radius: 4px;
                            cursor: pointer;
                            font-family: system-ui, -apple-system, sans-serif;
                            font-size: 16px;
                        }}
                        .select-btn:hover {{
                            background-color: #e4e4e4;
                        }}
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
                    wrappers.forEach(wrapper => {{
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
                    }});

                    // event listeners
                    document.getElementById('select-all-btn').addEventListener('click', () => {{
                        document.querySelectorAll('.sequence-checkbox').forEach(cb => {{
                            cb.checked = true;
                        }});
                    }});

                    document.getElementById('deselect-all-btn').addEventListener('click', () => {{
                        document.querySelectorAll('.sequence-checkbox').forEach(cb => {{
                            cb.checked = false;
                        }});
                    }});

                    document.getElementById('export-tsv-btn').addEventListener('click', exportSelectedAsTSV);
                }}

                // export selected data as TSV
                function exportSelectedAsTSV() {{
                    const selected = Array.from(document.querySelectorAll('.sequence-checkbox:checked'))
                        .map(cb => cb.getAttribute('data-contig'));
                    
                    if (selected.length === 0) {{
                        alert('Please select at least one sequence to export');
                        return;
                    }}
                    
                    // data for each selected contig
                    const exportData = [];
                    const headers = [
                        'Contig ID', // 'Length', 
                        'BLASTx Subject', 'BLASTx Organism', 'BLASTx Coverage', 'BLASTx Identity', 'BLASTx E-value',
                        'BLASTn Subject', 'BLASTn Coverage', 'BLASTn Identity', 'BLASTn E-value'
                    ];
                    
                    selected.forEach(contigId => {{
                        // viral hit data for this contig
                        const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${{contigId}}"]`);
                        
                        // extract information from the container
                        const infoPanel = wrapper.querySelector('.info-panel');
                        if (!infoPanel) return;
                        
                        // parse text content from the info panel
                        const infoText = infoPanel.textContent;
                        
                        // row object
                        const row = {{
                            'Contig ID': contigId
                        //'Length': extractValue(infoText, 'Length') || ''
                        }};
                        
                        // BLASTx data
                        if (infoText.includes('BLASTx Results')) {{
                            const blastxSection = infoText.split('BLASTx Results')[1].split('BLASTn Results')[0];
                            row['BLASTx Subject'] = extractValue(blastxSection, 'Subject') || '';
                            row['BLASTx Organism'] = extractValue(blastxSection, 'Organism') || '';
                            row['BLASTx Coverage'] = extractValue(blastxSection, 'Coverage') || '';
                            row['BLASTx Identity'] = extractValue(blastxSection, 'Identity') || '';
                            row['BLASTx E-value'] = extractValue(blastxSection, 'E-value') || '';
                        }}
                        
                        // BLASTn data
                        if (infoText.includes('BLASTn Results')) {{
                            const blastnSection = infoText.split('BLASTn Results')[1];
                            row['BLASTn Subject'] = extractValue(blastnSection, 'Subject') || '';
                            row['BLASTn Coverage'] = extractValue(blastnSection, 'Coverage') || '';
                            row['BLASTn Identity'] = extractValue(blastnSection, 'Identity') || '';
                            row['BLASTn E-value'] = extractValue(blastnSection, 'E-value') || '';
                        }}
                        
                        exportData.push(row);
                    }});
                    
                    // TSV content
                    let tsvContent = headers.join('\t') + '\n';
                    
                    exportData.forEach(row => {{
                        const rowValues = headers.map(header => row[header] || '');
                        tsvContent += rowValues.join('\t') + '\n';
                    }});
                    
                    // trigger download
                    const blob = new Blob([tsvContent], {{ type: 'text/tab-separated-values' }});
                    const url = URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'viralquest_export.tsv';
                    a.click();
                    URL.revokeObjectURL(url);
                }}

                // helper function to extract values from text content
                function extractValue(text, label) {{
                    const regex = new RegExp(label + ':\s*([^\n]+)');
                    const match = text.match(regex);
                    return match ? match[1].trim() : null;
                }}

                // run the initialization after visualizations are created
                const observer = new MutationObserver(function(mutations) {{
                    mutations.forEach(function(mutation) {{
                        if (mutation.type === 'childList' && 
                            document.querySelectorAll('.visualization-wrapper').length > 0 &&
                            document.getElementById('loading').style.display === 'none') {{
                            // visualizations are loaded, initialize our features
                            initializeSelectionAndExport();
                            // stop observing once initialization is complete
                            observer.disconnect();
                        }}
                    }});
                }});

                // start observing
                observer.observe(document.getElementById('visualization-container'), {{ childList: true, subtree: true }});

                // also try to initialize after a few seconds
                setTimeout(function() {{
                    if (document.querySelectorAll('.visualization-wrapper').length > 0) {{
                        initializeSelectionAndExport();
                        observer.disconnect();
                    }}
                }}, 3000);

                // initialize enhanced export functionality
                function initializeEnhancedExport() {{
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
                document.getElementById('select-all-btn').addEventListener('click', () => {{
                    document.querySelectorAll('.sequence-checkbox').forEach(cb => {{
                    cb.checked = true;
                    }});
                }});

                document.getElementById('deselect-all-btn').addEventListener('click', () => {{
                    document.querySelectorAll('.sequence-checkbox').forEach(cb => {{
                    cb.checked = false;
                    }});
                }});

                document.getElementById('export-tsv-btn').addEventListener('click', () => {{
                    const format = document.getElementById('export-format').value;
                    if (format === 'tsv') {{
                    exportSelectedAsTSV();
                    }} else if (format === 'svg') {{
                    exportSelectedAsGraphics('svg');
                    }} else if (format === 'pdf') {{
                    exportSelectedAsGraphics('pdf');
                    }}
                }});
                
                // customization dialog
                document.getElementById('customize-columns').addEventListener('click', showColumnsDialog);
                
                // close dialog when clicking overlay
                document.getElementById('dialog-overlay').addEventListener('click', hideColumnsDialog);
                }}

                // columns selection dialog
                function showColumnsDialog() {{
                    const columnsDialog = document.getElementById('columns-dialog');
                    const overlay = document.getElementById('dialog-overlay');
                
                // define all possible columns from jsonDB
                const possibleColumns = [
                    // viral Hits columns
                    {{ id: 'Sample_name', label: 'Sample Name', category: 'General', checked: true }},
                    {{ id: 'QueryID', label: 'Contig ID', category: 'General', checked: true }},
                    {{ id: 'BLASTx_Qlength', label: 'BLASTx Query Length', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTx_Slength', label: 'BLASTx Subj. Length', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTx_Cover', label: 'BLASTx Coverage', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTx_Ident', label: 'BLASTx Identity', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTx_evalue', label: 'BLASTx E-value', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTx_Subject_Title', label: 'BLASTx Subject Title', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTx_Organism_Name', label: 'BLASTx Organism', category: 'BLASTx', checked: true }},
                    {{ id: 'BLASTn_Qlength', label: 'BLASTn Query Length', category: 'BLASTn', checked: false }},
                    {{ id: 'BLASTn_Slength', label: 'BLASTn Subj. Length', category: 'BLASTn', checked: false }},
                    {{ id: 'BLASTn_Cover', label: 'BLASTn Coverage', category: 'BLASTn', checked: false }},
                    {{ id: 'BLASTn_Ident', label: 'BLASTn Identity', category: 'BLASTn', checked: false }},
                    {{ id: 'BLASTn_evalue', label: 'BLASTn E-value', category: 'BLASTn', checked: false }},
                    {{ id: 'BLASTn_Subject_Title', label: 'BLASTn Subject Title', category: 'BLASTn', checked: false }},
                    {{ id: 'TaxId', label: 'Taxonomy ID', category: 'Taxonomy', checked: false }},
                    {{ id: 'ScientificName', label: 'Scientific Name', category: 'Taxonomy', checked: true }},
                    {{ id: 'Clade', label: 'Clade', category: 'Taxonomy', checked: false }},
                    {{ id: 'Kingdom', label: 'Kingdom', category: 'Taxonomy', checked: true }},
                    {{ id: 'Phylum', label: 'Phylum', category: 'Taxonomy', checked: true }},
                    {{ id: 'Class', label: 'Class', category: 'Taxonomy', checked: true }},
                    {{ id: 'Order', label: 'Order', category: 'Taxonomy', checked: true }},
                    {{ id: 'Family', label: 'Family', category: 'Taxonomy', checked: true }},
                    {{ id: 'Genus', label: 'Genus', category: 'Taxonomy', checked: true }},
                    {{ id: 'Species', label: 'Species', category: 'Taxonomy', checked: true }},
                    {{ id: 'Genome', label: 'Genome Type', category: 'Taxonomy', checked: true }}
                ];
                
                // group columns by category
                const categories = [...new Set(possibleColumns.map(col => col.category))];
                
                // dialog content
                let dialogContent = `
                    <h3>Customize Export Columns</h3>
                    <p>Select the columns to include in your TSV export:</p>
                `;
                
                categories.forEach(category => {{
                    dialogContent += `<h4>${{category}}</h4>`;
                    dialogContent += '<div class="columns-grid">';
                    
                    possibleColumns.filter(col => col.category === category).forEach(col => {{
                        dialogContent += `
                            <div class="column-item">
                            <input type="checkbox" id="col-${{col.id}}" data-column="${{col.id}}" ${{col.checked ? 'checked' : ''}}>
                            <label for="col-${{col.id}}">${{col.label}}</label>
                            </div>
                        `;
                    }});
                    
                    dialogContent += '</div>';
                }});
                
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
                document.getElementById('select-all-columns').addEventListener('click', () => {{
                    columnsDialog.querySelectorAll('input[type="checkbox"]').forEach(cb => {{
                        cb.checked = true;
                    }});
                }});    
                
                document.getElementById('deselect-all-columns').addEventListener('click', () => {{
                    columnsDialog.querySelectorAll('input[type="checkbox"]').forEach(cb => {{
                        cb.checked = false;
                    }});
                }});
                
                document.getElementById('apply-columns').addEventListener('click', () => {{
                    // save selected columns to localStorage for persistence
                    const selectedColumns = Array.from(columnsDialog.querySelectorAll('input[type="checkbox"]'))
                        .filter(cb => cb.checked)
                        .map(cb => cb.dataset.column);
                    
                    localStorage.setItem('vqSelectedColumns', JSON.stringify(selectedColumns));
                    hideColumnsDialog();
                }});
                
                document.getElementById('cancel-columns').addEventListener('click', hideColumnsDialog);
                
                // dialog and overlay
                columnsDialog.style.display = 'block';
                overlay.style.display = 'block';
                }}

                // hide columns selection dialog
                function hideColumnsDialog() {{
                    document.getElementById('columns-dialog').style.display = 'none';
                    document.getElementById('dialog-overlay').style.display = 'none';
                }}

                // function to export selected contigs as TSV with data from jsonDB
                function exportSelectedAsTSV() {{
                    const selected = Array.from(document.querySelectorAll('.sequence-checkbox:checked'))
                        .map(cb => cb.getAttribute('data-contig'));
                    
                    if (selected.length === 0) {{
                        alert('Please select at least one sequence to export');
                        return;
                    }}
                    
                    // selected columns from localStorage or use defaults
                    let selectedColumns;
                    try {{
                        selectedColumns = JSON.parse(localStorage.getItem('vqSelectedColumns')) || [
                        'QueryID', 'Sample_name', 'BLASTx_Qlength', 'BLASTx_Cover', 'BLASTx_Ident', 
                        'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name',
                        'ScientificName', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
                        ];
                    }} catch (e) {{
                        // default columns if parsing fails
                        selectedColumns = [
                        'QueryID', 'Sample_name', 'BLASTx_Qlength', 'BLASTx_Cover', 'BLASTx_Ident', 
                        'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name',
                        'ScientificName', 'Kingdom', 'Family', 'Genus', 'Species'
                        ];
                    }}
                    
                    // column labels (for header row)
                    const columnLabels = {{
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
                    }};
                    
                    // header row
                    const headers = selectedColumns.map(col => columnLabels[col] || col);
                    
                    // data for selected contigs directly from jsonDB
                    const rows = [];
                    
                    selected.forEach(contigId => {{
                        // find contigs Viral_Hits data
                        const hit = jsonData.Viral_Hits.find(h => h.QueryID === contigId);
                        if (hit) {{
                        const row = {{}};
                        
                        // values for selected columns
                        selectedColumns.forEach(col => {{
                            row[columnLabels[col] || col] = hit[col] !== undefined ? hit[col] : '';
                        }});
                        
                        rows.push(row);
                        }}
                    }});
                    
                    // TSV content
                    let tsvContent = headers.join('\t') + '\n';
                    
                    rows.forEach(row => {{
                        const rowValues = headers.map(header => row[header] || '');
                        tsvContent += rowValues.join('\t') + '\n';
                    }});
                    
                    // trigger download
                    const blob = new Blob([tsvContent], {{ type: 'text/tab-separated-values' }});
                    const url = URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'viralquest_export.tsv';
                    a.click();
                    URL.revokeObjectURL(url);
                }}

                // export selected visualizations as SVG or PDF
                function exportSelectedAsGraphics(format) {{
                    const selected = Array.from(document.querySelectorAll('.sequence-checkbox:checked'))
                        .map(cb => cb.getAttribute('data-contig'));
                    
                    if (selected.length === 0) {{
                        alert('Please select at least one visualization to export');
                        return;
                    }}
                    
                    if (format === 'svg') {{
                        // Export as individual SVG files or as a combined SVG
                        if (selected.length === 1) {{
                        // Single SVG export
                        exportSingleSVG(selected[0]);
                        }} else {{
                        // multiple SVGs - ask if user wants individual files or a combined file
                        if (confirm('Export as individual SVG files? Click Cancel for a combined file.')) {{
                            selected.forEach(contigId => exportSingleSVG(contigId));
                        }} else {{
                            exportCombinedSVG(selected);
                        }}
                        }}
                    }} else if (format === 'pdf') {{
                        // for PDF we need to use a library like jsPDF
                        // check if jsPDF is loaded
                        if (typeof jspdf === 'undefined') {{
                        // load jsPDF dynamically
                        const script = document.createElement('script');
                        script.src = 'https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js';
                        script.onload = function() {{
                            // also load svg2pdf
                            const svg2pdf = document.createElement('script');
                            svg2pdf.src = 'https://cdnjs.cloudflare.com/ajax/libs/svg2pdf.js/2.2.1/svg2pdf.min.js';
                            svg2pdf.onload = function() {{
                            exportAsPDF(selected);
                            }};
                            document.head.appendChild(svg2pdf);
                        }};
                        document.head.appendChild(script);
                        }} else {{
                        // jsPDF is already loaded
                        exportAsPDF(selected);
                        }}
                    }}
                }}

                // wxport a single visualization as SVG
                function exportSingleSVG(contigId) {{
                    const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${{contigId}}"]`);
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
                    const blob = new Blob([svgData], {{ type: 'image/svg+xml' }});
                    const url = URL.createObjectURL(blob);
                    
                    // trigger download
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = `${{contigId}}.svg`;
                    a.click();
                    URL.revokeObjectURL(url);
                }}

                // export multiple visualizations as a combined SVG
                function exportCombinedSVG(contigIds) {{
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
                    contigIds.forEach((contigId, index) => {{
                        const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${{contigId}}"]`);
                        if (!wrapper) return;
                        
                        const svg = wrapper.querySelector('svg');
                        if (!svg) return;
                        
                        // visualization content 
                        const mainGroup = svg.querySelector('g');
                        if (!mainGroup) return;
                        
                        // clone it
                        const groupClone = mainGroup.cloneNode(true);
                        
                        // position the group at the correct vertical position
                        groupClone.setAttribute('transform', `translate(20, ${{totalHeight}})`);
                        
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
                    }});
                    
                    // set dimensions of the combined SVG
                    combinedSvg.setAttribute('width', maxWidth);
                    combinedSvg.setAttribute('height', totalHeight);
                    combinedSvg.setAttribute('viewBox', `0 0 ${{maxWidth}} ${{totalHeight}}`);
                    
                    // styling
                    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
                    style.textContent = `
                        text {{ font-family: system-ui, -apple-system, sans-serif; }}
                        .genome-line {{ stroke: #2c3e50; stroke-width: 2; }}
                        .scale-axis path, .scale-axis line {{ stroke: #7f8c8d; }}
                    `;
                    combinedSvg.insertBefore(style, combinedSvg.firstChild);
                    
                    // SVG content
                    const svgData = new XMLSerializer().serializeToString(combinedSvg);
                    
                    // downloadable blob
                    const blob = new Blob([svgData], {{ type: 'image/svg+xml' }});
                    const url = URL.createObjectURL(blob);
                    
                    // trigger download
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'viralquest_combined.svg';
                    a.click();
                    URL.revokeObjectURL(url);
                }}

                // enhance SVG before export
                function enhanceSVGForExport(svgElement) {{
                    // CSS styles directly to the SVG
                    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
                    style.textContent = `
                        text {{ font-family: system-ui, -apple-system, sans-serif; }}
                        .genome-line {{ stroke: #2c3e50; stroke-width: 2; }}
                        .scale-axis path, .scale-axis line {{ stroke: #7f8c8d; }}
                        .pfam-domain {{ stroke: #333; stroke-width: 0.5; }}
                    `;
                    svgElement.insertBefore(style, svgElement.firstChild);
                    
                    // title and description
                    const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
                    title.textContent = 'ViralQuest Visualization';
                    svgElement.insertBefore(title, svgElement.firstChild);
                    
                    // ensure elements have proper styling attributes
                    const paths = svgElement.querySelectorAll('path');
                    paths.forEach(path => {{
                        if (!path.getAttribute('stroke') && !path.getAttribute('style').includes('stroke')) {{
                        path.setAttribute('stroke', '#333');
                        path.setAttribute('stroke-width', '0.5');
                        }}
                    }});
                    
                    return svgElement;
                }}

                // Function to export as PDF
                // function exportAsPDF(contigIds) {{
                //     // Check if required libraries are loaded                
                //     if (typeof jspdf === 'undefined' || typeof svg2pdf === 'undefined') {{
                //         alert('Required libraries for PDF export are not available. Please try again.');
                //         return;
                //     }}
                    
                //     // Create a new PDF document
                //     const {{ jsPDF }} = jspdf;
                //     const doc = new jsPDF('landscape', 'pt');
                //     let currentPage = 1;
                    
                //     // Process each selected contig
                //     const processNextContig = (index) => {{
                //         if (index >= contigIds.length) {{
                //         // All contigs processed, save the PDF
                //         doc.save('viralquest_export.pdf');
                //         return;
                //         }}
                        
                //         const contigId = contigIds[index];
                //         const wrapper = document.querySelector(`.visualization-wrapper[data-contig="${{contigId}}"]`);
                //         if (!wrapper) {{
                //         // Skip this one and move to next
                //         processNextContig(index + 1);
                //         return;
                //         }}
                        
                //         const svg = wrapper.querySelector('svg');
                //         if (!svg) {{
                //         // Skip this one and move to next
                //         processNextContig(index + 1);
                //         return;
                //         }}
                        
                //         // Clone the SVG for export
                //         const svgClone = svg.cloneNode(true);
                //         enhanceSVGForExport(svgClone);
                        
                //         // Add a new page if not the first contig
                //         if (index > 0) {{
                //         doc.addPage();
                //         currentPage++;
                //         }}
                        
                //         // Add title to the page
                //         doc.setFontSize(16);
                //         doc.text(`Contig: ${{contigId}}`, 40, 40);
                        
                //         // Convert SVG to PDF
                //         const svgElement = svgClone;
                //         svgElement.setAttribute('width', doc.internal.pageSize.getWidth() - 80);
                //         svgElement.setAttribute('height', doc.internal.pageSize.getHeight() - 100);
                        
                //         // Use svg2pdf to add the SVG to the PDF
                //         svg2pdf(doc, svgElement, {{
                //         xOffset: 40,
                //         yOffset: 60,
                //         width: doc.internal.pageSize.getWidth() - 80,
                //         height: doc.internal.pageSize.getHeight() - 100
                //         }})
                //         .then(() => {{
                //         // Add page number
                //         doc.setFontSize(10);
                //         doc.text(`Page ${{currentPage}} of ${{contigIds.length}}`, 
                //                 doc.internal.pageSize.getWidth() - 100, 
                //                 doc.internal.pageSize.getHeight() - 30);
                        
                //         // Process next contig
                //         processNextContig(index + 1);
                //         }});
                //     }};
                    
                //     // Start processing with the first contig
                //     processNextContig(0);
                // }}            
                // <option value="pdf">PDF (Graphics)</option>


        
                function enhanceExistingExport() {{
                    // check if basic export panel exists
                    const existingPanel = document.querySelector('.export-panel');
                    if (existingPanel) {{
                        // replace with enhanced version
                        initializeEnhancedExport();
                    }} else {{
                        setTimeout(enhanceExistingExport, 500);
                    }}
                }}

                // call function
                document.addEventListener('DOMContentLoaded', () => {{
                    // check basic export panel
                    setTimeout(enhanceExistingExport, 1000);
                }});

        </script>
    </body>
    </html>
    """

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)


