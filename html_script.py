import csv
import json

def csv_to_json(csv_file_path1, csv_file_path2):
    with open(csv_file_path1, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        viral_hits = [row for row in csv_reader]
    
    with open(csv_file_path2, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        hmm_hits = [row for row in csv_reader]

    return {"Viral_Hits": viral_hits, "HMM_hits": hmm_hits}

#json_result = csv_to_json('break_allHits.csv', 'break_hmm.csv')
#print(json.dumps(json_result, indent=4))

def generate_html_report(json_data, output_file="report.html"):
    viral_hits_json = json.dumps(json_data["Viral_Hits"], indent=4)
    hmm_hits_json = json.dumps(json_data["HMM_hits"], indent=4)

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>ViralQuest Report</title>
        <script src="https://d3js.org/d3.v7.min.js"></script>
        <style>
            :root {{
                --primary-color: #0a3f6b;
                --secondary-color: #69b3a2;
                --text-color: #333;
                --background-color: #f5f5f5;
                --sidebar-width: 220px;
            }}
            
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                color: var(--text-color);
                background-color: var(--background-color);
            }}
            
            .sidenav {{
                height: 100vh;
                width: var(--sidebar-width);
                position: fixed;
                z-index: 1;
                top: 0;
                left: 0;
                background-color: var(--primary-color);
                overflow-x: hidden;
                padding-top: 40px;
                box-shadow: 2px 0 5px rgba(0,0,0,0.1);
            }}
            
            .sidenav a {{
                padding: 15px 25px;
                text-decoration: none;
                font-size: 16px;
                color: #fff;
                display: block;
                transition: all 0.3s ease;
                border-left: 4px solid transparent;
            }}
            
            .sidenav a:hover {{
                background-color: rgba(255,255,255,0.1);
                border-left: 4px solid var(--secondary-color);
            }}
            
            .main {{
                margin-left: var(--sidebar-width);
                padding: 40px;
                min-height: 100vh;
            }}
            
            h1 {{
                color: var(--primary-color);
                margin-bottom: 30px;
                font-size: 2.5em;
            }}
            
            h2 {{
                color: var(--primary-color);
                margin-bottom: 20px;
                font-size: 1.8em;
            }}
            
            section {{
                background: white;
                padding: 30px;
                margin: 30px 0;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            
            #my_dataviz {{
                width: 100%;
                height: 500px;
                margin-top: 20px;
            }}
            
            .button-group {{
                display: flex;
                gap: 10px;
                margin-bottom: 20px;
            }}
            
            button {{
                padding: 12px 24px;
                background-color: var(--primary-color);
                color: white;
                border: none;
                border-radius: 6px;
                cursor: pointer;
                font-size: 14px;
                transition: all 0.3s ease;
            }}
            
            button:hover {{
                background-color: #0d4d82;
                transform: translateY(-1px);
            }}
            
            .chart-container {{
                position: relative;
                width: 100%;
                height: 100%;
            }}

            /* Table Styles */
            .data-table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
                font-size: 14px;
            }}

            .data-table th {{
                background-color: var(--primary-color);
                color: white;
                padding: 12px;
                text-align: left;
                position: sticky;
                top: 0;
            }}

            .data-table td {{
                padding: 10px;
                border-bottom: 1px solid #ddd;
            }}

            .data-table tr:nth-child(even) {{
                background-color: #f8f9fa;
            }}

            .data-table tr:hover {{
                background-color: #f0f0f0;
            }}

            .table-container {{
                max-height: 500px;
                overflow-y: auto;
                margin-top: 20px;
            }}

            .table-controls {{
                margin-bottom: 20px;
            }}

            .search-input {{
                padding: 8px;
                width: 200px;
                border: 1px solid #ddd;
                border-radius: 4px;
                margin-right: 10px;
            }}
        </style>
    </head>
    <body>
        <div class="sidenav">
            <a href="#viral-distribution">Distribution Analysis</a>
            <a href="#viral-hits">Viral Hits</a>
            <a href="#hmm-hits">HMM Hits</a>
        </div>

        <div class="main">
            <h1>ViralQuest Analysis Report</h1>
            
            <section id="viral-distribution">
                <h2>Viral Contigs Family Distribution</h2>
                <div class="button-group">
                    <button onclick="update(dataViralHits)">Show Viral Hits</button>
                    <button onclick="update(dataHMMHits)">Show HMM Hits</button>
                </div>
                <div class="chart-container">
                    <div id="my_dataviz"></div>
                </div>
            </section>

            <section id="viral-hits">
                <h2>Filtered Viral Hits</h2>
                <div class="table-controls">
                    <input type="text" id="viral-search" class="search-input" placeholder="Search viral hits...">
                </div>
                <div class="table-container" id="viral-hits-table"></div>
            </section>

            <section id="hmm-hits">
                <h2>HMM Hits</h2>
                <div class="table-controls">
                    <input type="text" id="hmm-search" class="search-input" placeholder="Search HMM hits...">
                </div>
                <div class="table-container" id="hmm-hits-table"></div>
            </section>
        </div>

        <script>
            // Embed the data
            const viralHits = {viral_hits_json};
            const hmmHits = {hmm_hits_json};

            // Prepare data for visualization
            function prepareChartData(data, familyKey) {{
                return Object.entries(data.reduce((acc, row) => {{
                    const family = row[familyKey] || "Unknown";
                    acc[family] = (acc[family] || 0) + 1;
                    return acc;
                }}, {{}}))
                .map(([group, value]) => {{ return {{ group, value }}; }});
            }}

            const dataViralHits = prepareChartData(viralHits, 'BLASTx_Family');
            const dataHMMHits = prepareChartData(hmmHits, 'ViralFamily');

            // Function to calculate dynamic font size
            function calculateFontSize(dataLength) {{
                const maxFontSize = 14;
                const minFontSize = 8;
                const maxItems = 20;
                
                if (dataLength <= 5) {{
                    return maxFontSize;
                }} else if (dataLength >= maxItems) {{
                    return minFontSize;
                }} else {{
                    return maxFontSize - ((dataLength - 5) * (maxFontSize - minFontSize) / (maxItems - 5));
                }}
            }}

            // Set up the chart dimensions
            const margin = {{ top: 40, right: 30, bottom: 100, left: 60 }};
            const container = document.getElementById("my_dataviz");
            const width = container.clientWidth - margin.left - margin.right;
            const height = container.clientHeight - margin.top - margin.bottom;

            // Create SVG
            const svg = d3.select("#my_dataviz")
                .append("svg")
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("viewBox", `0 0 ${{container.clientWidth}} ${{container.clientHeight}}`)
                .append("g")
                .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

            // Initialize scales
            const x = d3.scaleBand().range([0, width]).padding(0.2);
            const y = d3.scaleLinear().range([height, 0]);

            // Add axes
            const xAxis = svg.append("g")
                .attr("transform", `translate(0,${{height}})`);
            
            const yAxis = svg.append("g")
                .attr("class", "myYaxis");

            // Add Y-axis label
            svg.append("text")
                .attr("class", "y-axis-label")
                .attr("transform", "rotate(-90)")
                .attr("y", 0 - margin.left)
                .attr("x", 0 - (height / 2))
                .attr("dy", "1em")
                .style("text-anchor", "middle")
                .text("Count");

            // Update function with dynamic font sizing
            function update(data) {{
                // Calculate font size based on data length
                const fontSize = calculateFontSize(data.length);
                
                // Update scales
                x.domain(data.map(d => d.group));
                y.domain([0, d3.max(data, d => d.value)]);

                // Update X axis with dynamic font size
                xAxis.transition()
                    .duration(1000)
                    .call(d3.axisBottom(x))
                    .selectAll("text")
                    .style("text-anchor", "end")
                    .style("font-size", `${{fontSize}}px`)
                    .attr("dx", "-.8em")
                    .attr("dy", ".15em")
                    .attr("transform", "rotate(-45)");

                // Update Y axis
                yAxis.transition()
                    .duration(1000)
                    .call(d3.axisLeft(y));

                // Update bars
                const bars = svg.selectAll("rect")
                    .data(data);

                // Remove old bars
                bars.exit().remove();

                // Add new bars with transition
                bars.enter()
                    .append("rect")
                    .merge(bars)
                    .transition()
                    .duration(1000)
                    .attr("x", d => x(d.group))
                    .attr("y", d => y(d.value))
                    .attr("width", x.bandwidth())
                    .attr("height", d => height - y(d.value))
                    .attr("fill", "#69b3a2");
            }}

            // Handle window resize
            window.addEventListener('resize', () => {{
                const container = document.getElementById("my_dataviz");
                svg.attr("viewBox", `0 0 ${{container.clientWidth}} ${{container.clientHeight}}`);
                
                // Update scales
                x.range([0, container.clientWidth - margin.left - margin.right]);
                y.range([container.clientHeight - margin.top - margin.bottom, 0]);
                
                // Update current data
                update(dataViralHits);
            }});

            // Initialize with Viral Hits data
            update(dataViralHits);

            // Table creation functions
            function createTable(data, containerId) {{
                if (data.length === 0) return;
                
                const container = document.getElementById(containerId);
                const keys = Object.keys(data[0]);
                
                const table = document.createElement('table');
                table.className = 'data-table';
                
                // Create header
                const thead = document.createElement('thead');
                const headerRow = document.createElement('tr');
                keys.forEach(key => {{
                    const th = document.createElement('th');
                    th.textContent = key;
                    headerRow.appendChild(th);
                }});
                thead.appendChild(headerRow);
                table.appendChild(thead);
                
                // Create body
                const tbody = document.createElement('tbody');
                data.forEach(row => {{
                    const tr = document.createElement('tr');
                    keys.forEach(key => {{
                        const td = document.createElement('td');
                        td.textContent = row[key] || '';
                        tr.appendChild(td);
                    }});
                    tbody.appendChild(tr);
                }});
                table.appendChild(tbody);
                
                container.innerHTML = '';
                container.appendChild(table);
            }}

            // Table filter function
            function filterTable(data, searchTerm) {{
                if (!searchTerm) return data;
                searchTerm = searchTerm.toLowerCase();
                return data.filter(row => 
                    Object.values(row).some(value => 
                        String(value).toLowerCase().includes(searchTerm)
                    )
                );
            }}

            // Initialize tables
            createTable(viralHits, 'viral-hits-table');
            createTable(hmmHits, 'hmm-hits-table');

            // Add search functionality
            document.getElementById('viral-search').addEventListener('input', (e) => {{
                const filtered = filterTable(viralHits, e.target.value);
                createTable(filtered, 'viral-hits-table');
            }});

            document.getElementById('hmm-search').addEventListener('input', (e) => {{
                const filtered = filterTable(hmmHits, e.target.value);
                createTable(filtered, 'hmm-hits-table');
            }});
        </script>
    </body>
    </html>
    """


    with open(output_file, "w", encoding="utf-8") as file:
        file.write(html_content)



# Example Usage
json_data = csv_to_json("break_allHits.csv", "break_hmm.csv")
generate_html_report(json_data, output_file="report.html")
print("HTML report generated successfully.")