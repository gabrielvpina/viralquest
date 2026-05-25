import os, glob, re, time
import pandas as pd

from .data import viralInfo_table

def remove_thinking_tags(text):

    if hasattr(text, 'content'):
        text_content = text.content
    elif isinstance(text, dict) and 'content' in text:
        text_content = text['content']
    else:
        text_content = str(text)
    
    # remove <thinking>...</thinking> blocks
    cleaned_text = re.sub(r'<thinking>.*?</thinking>', '', text_content, flags=re.DOTALL)
    cleaned_text = re.sub(r'<think>.*?</think>', '', cleaned_text, flags=re.DOTALL)
    # also handle variations 
    cleaned_text = re.sub(r'<[Tt]hinking>.*?</[Tt]hinking>', '', cleaned_text, flags=re.DOTALL)
    cleaned_text = re.sub(r'<[Tt]hink>.*?</[Tt]hink>', '', cleaned_text, flags=re.DOTALL)
    # remove any extra whitespace that might be left
    cleaned_text = re.sub(r'\s+', ' ', cleaned_text).strip()
    return cleaned_text



_last_api_call_time = 0

def get_llm_model(args):

    """
    Initialize and return the appropriate LLM based on user arguments.
    For 'ollama' type, this function will return None as it's handled directly.
    """
    if args.model_type == "ollama":
        # no langchain llm object needed for direct ollama usage
        return None
    elif args.model_type == "openai":
        from langchain_openai import ChatOpenAI
        if not args.api_key:
            raise ValueError("API key is required for OpenAI models")
        os.environ["OPENAI_API_KEY"] = args.api_key
        return ChatOpenAI(model_name=args.model_name, temperature=0.3)
    elif args.model_type == "anthropic":
        from langchain_anthropic import ChatAnthropic
        if not args.api_key:
            raise ValueError("API key is required for Anthropic models")
        os.environ["ANTHROPIC_API_KEY"] = args.api_key
        return ChatAnthropic(model_name=args.model_name, temperature=0.3)
    elif args.model_type == "google":
        from langchain_google_genai import ChatGoogleGenerativeAI
        if not args.api_key:
            raise ValueError("API key is required for Google models")
        os.environ["GOOGLE_API_KEY"] = args.api_key
        return ChatGoogleGenerativeAI(model=args.model_name, temperature=0.3)
    else:
        raise ValueError(f"Unsupported model type: {args.model_type}")

def analyze_viral_sequences(args):
    
    from langchain_core.prompts import ChatPromptTemplate
    import ollama
    import json
    import re

    resultFolder = os.path.normpath(args.outdir)
    fileJSON = os.path.join(resultFolder, "*.json")
    viralSeqs_path = glob.glob(fileJSON)
    json_file_path = viralSeqs_path[0]

    output_json_path = json_file_path
    
    # json data
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    # viral taxonomy data
    taxonomy_df = viralInfo_table
    
    # start llm - initialize langchain model only if not ollama
    langchain_llm_model = None
    if args.model_type != "ollama":
        try:
            langchain_llm_model = get_llm_model(args)
            if langchain_llm_model is None:
                raise ValueError(f"Model could not be initialized for type: {args.model_type}")
            # print(f"Successfully loaded {args.model_type} model: {args.model_name}")
        except Exception as e:
            print(f"Error initializing Langchain model: {e}")
            return None
    elif args.model_type == "ollama":
        try:
            # check if ollama server is accessible and model exists
            ollama.show(args.model_name)
            # print(f"Successfully connected to Ollama. Using model: {args.model_name}")
        except Exception as e:
            print(f"Error connecting to Ollama or model '{args.model_name}' not found/available: {e}")
            print(f"Please ensure Ollama is running and the model is downloaded (e.g., 'ollama pull {args.model_name}').")
            return None
    
    # template for langchain models (openai, anthropic, google)
    langchain_template = """
        You are a virology expert analyzing potential viral sequences from bioinformatics data. Be skeptical in your analysis.

        ## Input Data:
        - BLASTx and BLASTn hits with organism taxonomy
        - HMMER hits (RVDB, Vfam, eggNOG, Pfam)
        - Viral taxonomy information (order/family/genus)

        Sequence data: 
        {sequence_data}

        Taxonomy info: 
        {taxonomy_info}

        ## Required Outputs:

        ### 1. Text Analysis (max 200 words)
        Analyze taxonomy, host range, geographic circulation, and classify as:
        - Known virus (>90% identity + >70% coverage)
        - Novel virus 
        - Non-viral sequence
        Include demarcation criteria for the taxonomic group.
        
        ### 2. VQ_SCORE (0-100)
        Calculate based on weighted criteria:

        **BLASTn identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts  
        - <70% identity or <50% coverage = 8-15 pts
        - No data/non-viral hit = 0-8 pts

        **BLASTx identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts
        - <70% identity or <50% coverage = 0-15 pts

        **HMM domain detection (40%)**
        - Strong viral domains (score >100) = 30-40 pts
        - Moderate viral domains = 20-30 pts
        - Weak/no domains = 0-20 pts

        **Score interpretation:**
        - 90-100: Known viral sequences
        - 70-90: High confidence viral
        - 50-70: Likely novel viruses  
        - 30-50: Low information
        - 0-30: Low viral probability
        
        ## Format Requirements:
        1. Maximum 200 words for the text analysis
        2. Plain text only for analysis
        3. No introduction, greeting, or personal commentary
        4. No reference to prompt or question
        5. Avoid using Markdown header elements in the text
        8. Don't provide any reference about VQ SCORE calculation in the text
        9. Handle missing data gracefully

        ## Example Output:
        [TEXT_ANALYSIS]
        The sequence shows 95% amino acid identity to Hepatitis B virus polymerase with 85% coverage
        via BLASTx. HMMER detected strong viral replication domains (score 150). This represents a
        known HBV strain with typical host range in humans and circulation in Asia-Pacific regions.
        Hepadnaviridae demarcation requires >80% nucleotide identity for species classification.
        [/TEXT_ANALYSIS]

        [VQ_SCORE]
        92
        [/VQ_SCORE]

        ## CRITICAL FORMATTING REQUIREMENTS:
        - The VQ_SCORE must be a plain integer between 0-100 with NO quotes, NO decimals, NO additional text
        - WRONG formats: [VQ_SCORE]"85"[/VQ_SCORE] or [VQ_SCORE]85.0[/VQ_SCORE] or [VQ_SCORE]about 85[/VQ_SCORE]
    """

    # prompts for direct ollama usage
    ollama_system_prompt = """
        You are a virology expert analyzing potential viral sequences from bioinformatics data. Be skeptical in your analysis.

        ## Input Data:
        - BLASTx and BLASTn hits with organism taxonomy
        - HMMER hits (RVDB, Vfam, eggNOG, Pfam)
        - Viral taxonomy information (order/family/genus)

        ## Required Outputs:

        ### 1. Text Analysis (max 200 words)
        Analyze taxonomy, host range, geographic circulation, and classify as:
        - Known virus (>90% identity + >70% coverage)
        - Novel virus 
        - Non-viral sequence
        Include demarcation criteria for the taxonomic group.
        
        ### 2. VQ_SCORE (0-100)
        Calculate based on weighted criteria:

        **BLASTn identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts  
        - <70% identity or <50% coverage = 8-15 pts
        - No data/non-viral hit = 0-8 pts

        **BLASTx identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts
        - <70% identity or <50% coverage = 0-15 pts

        **HMM domain detection (40%)**
        - Strong viral domains (score >100) = 30-40 pts
        - Moderate viral domains = 20-30 pts
        - Weak/no domains = 0-20 pts

        **Score interpretation:**
        - 90-100: Known viral sequences
        - 70-90: High confidence viral
        - 50-70: Likely novel viruses  
        - 30-50: Low information
        - 0-30: Low viral probability
        
        ## Format Requirements:
        1. Maximum 200 words for the text analysis
        2. Plain text only for analysis
        3. No introduction, greeting, or personal commentary
        4. No reference to prompt or question
        5. Avoid using Markdown header elements in the text
        8. Don't provide any reference about VQ SCORE calculation in the text
        9. Handle missing data gracefully

        ## Example Output:
        [TEXT_ANALYSIS]
        The sequence shows 95% amino acid identity to Hepatitis B virus polymerase with 85% coverage
        via BLASTx. HMMER detected strong viral replication domains (score 150). This represents a
        known HBV strain with typical host range in humans and circulation in Asia-Pacific regions.
        Hepadnaviridae demarcation requires >80% nucleotide identity for species classification.
        [/TEXT_ANALYSIS]

        [VQ_SCORE]
        92
        [/VQ_SCORE]

        ## CRITICAL FORMATTING REQUIREMENTS:
        - The VQ_SCORE must be a plain integer between 0-100 with NO quotes, NO decimals, NO additional text
        - WRONG formats: [VQ_SCORE]"85"[/VQ_SCORE] or [VQ_SCORE]85.0[/VQ_SCORE] or [VQ_SCORE]about 85[/VQ_SCORE]
        - Don't write any reference of [VQ_SCORE] in [TEXT_ANALYSIS] section.

    """

    ollama_user_prompt_template = """
        Sequence data:
        {sequence_data}

        Taxonomy info:
        {taxonomy_info}
    """




    ################################### NO TAXONOMY

    langchain_template_noTax = """
        You are a virology expert analyzing potential viral sequences from bioinformatics data. Be skeptical in your analysis.

        ## Input Data:
        - BLASTx and BLASTn hits with organism taxonomy
        - HMMER hits (RVDB, Vfam, eggNOG, Pfam)

        Sequence data: 
        {sequence_data}

        ## Required Outputs:

        ### 1. Text Analysis (max 200 words)
        Analyze taxonomy, host range, geographic circulation, and classify as:
        - Known virus (>90% identity + >70% coverage)
        - Novel virus 
        - Non-viral sequence
        Include demarcation criteria for the taxonomic group.
        
        ### 2. VQ_SCORE (0-100)
        Calculate based on weighted criteria:

        **BLASTn identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts  
        - <70% identity or <50% coverage = 8-15 pts
        - No data/non-viral hit = 0-8 pts

        **BLASTx identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts
        - <70% identity or <50% coverage = 0-15 pts

        **HMM domain detection (40%)**
        - Strong viral domains (score >100) = 30-40 pts
        - Moderate viral domains = 20-30 pts
        - Weak/no domains = 0-20 pts

        **Score interpretation:**
        - 90-100: Known viral sequences
        - 70-90: High confidence viral
        - 50-70: Likely novel viruses  
        - 30-50: Low information
        - 0-30: Low viral probability
        
        ## Format Requirements:
        1. Maximum 200 words for the text analysis
        2. Plain text only for analysis
        3. No introduction, greeting, or personal commentary
        4. No reference to prompt or question
        5. Avoid using Markdown header elements in the text
        8. Don't provide any reference about VQ SCORE calculation in the text
        9. Handle missing data gracefully

        ## Example Output:
        [TEXT_ANALYSIS]
        The sequence shows 95% amino acid identity to Hepatitis B virus polymerase with 85% coverage
        via BLASTx. HMMER detected strong viral replication domains (score 150). These domains and 
        BLASTn and BLASTx results indicate a problable viral sequence.
        [/TEXT_ANALYSIS]

        [VQ_SCORE]
        92
        [/VQ_SCORE]

        ## CRITICAL FORMATTING REQUIREMENTS:
        - The VQ_SCORE must be a plain integer between 0-100 with NO quotes, NO decimals, NO additional text
        - WRONG formats: [VQ_SCORE]"85"[/VQ_SCORE] or [VQ_SCORE]85.0[/VQ_SCORE] or [VQ_SCORE]about 85[/VQ_SCORE]
    """

    # prompts for direct ollama usage
    ollama_system_prompt_noTax = """

        You are a virology expert analyzing potential viral sequences from bioinformatics data. Be skeptical in your analysis.

        ## Input Data:
        - BLASTx and BLASTn hits with organism taxonomy
        - HMMER hits (RVDB, Vfam, eggNOG, Pfam)

        ## Required Outputs:

        ### 1. Text Analysis (max 200 words)
        Analyze taxonomy, host range, geographic circulation, and classify as:
        - Known virus (>90% identity + >70% coverage)
        - Novel virus 
        - Non-viral sequence
        Include demarcation criteria for the taxonomic group.
        
        ### 2. VQ_SCORE (0-100)
        Calculate based on weighted criteria:

        **BLASTn identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts  
        - <70% identity or <50% coverage = 8-15 pts
        - No data/non-viral hit = 0-8 pts

        **BLASTx identity/coverage (30%)**
        - >90% identity + >70% coverage = 25-30 pts
        - 70-90% identity + >50% coverage = 15-25 pts
        - <70% identity or <50% coverage = 0-15 pts

        **HMM domain detection (40%)**
        - Strong viral domains (score >100) = 30-40 pts
        - Moderate viral domains = 20-30 pts
        - Weak/no domains = 0-20 pts

        **Score interpretation:**
        - 90-100: Known viral sequences
        - 70-90: High confidence viral
        - 50-70: Likely novel viruses  
        - 30-50: Low information
        - 0-30: Low viral probability
        
        ## Format Requirements:
        1. Maximum 200 words for the text analysis
        2. Plain text only for analysis
        3. No introduction, greeting, or personal commentary
        4. No reference to prompt or question
        5. Avoid using Markdown header elements in the text
        8. Don't provide any reference about VQ SCORE calculation in the text
        9. Handle missing data gracefully

        ## Example Output:
        [TEXT_ANALYSIS]
        The sequence shows 95% amino acid identity to Hepatitis B virus polymerase with 85% coverage
        via BLASTx. HMMER detected strong viral replication domains (score 150). These domains and 
        BLASTn and BLASTx results indicate a problable viral sequence.
        [/TEXT_ANALYSIS]

        [VQ_SCORE]
        92
        [/VQ_SCORE]

        ## CRITICAL FORMATTING REQUIREMENTS:
        - The VQ_SCORE must be a plain integer between 0-100 with NO quotes, NO decimals, NO additional text
        - WRONG formats: [VQ_SCORE]"85"[/VQ_SCORE] or [VQ_SCORE]85.0[/VQ_SCORE] or [VQ_SCORE]about 85[/VQ_SCORE]


    """

    ollama_user_prompt_template_noTax = """
    Input Data:
    {sequence_data}
    """


    def apply_api_rate_limit():
        """apply 40-second rate limiting for api calls"""
        _last_api_call_time = 0
        #global _last_api_call_time
        current_time = time.time()
        time_since_last_call = current_time - _last_api_call_time
        
        # if less than 20 seconds have passed, wait
        if time_since_last_call < 40:
            wait_time = 40 - time_since_last_call
            time.sleep(wait_time)
        
        # update the last call time
        _last_api_call_time = time.time()
    
    def parse_llm_response(response_text):
        """
        parse the llm response to extract text analysis and vq_score
        
        Args:
            response_text (str): the raw response from the llm
            
        Returns:
            tuple: (text_analysis, vq_score) where vq_score is int or None
        """
        text_analysis = ""
        vq_score = None
        
        try:
            # extract text analysis
            text_match = re.search(r'\[TEXT_ANALYSIS\](.*?)\[/TEXT_ANALYSIS\]', response_text, re.DOTALL)
            if text_match:
                text_analysis = text_match.group(1).strip()
            else:
                # fallback: if tags not found, try to extract everything before VQ_SCORE
                score_match = re.search(r'\[VQ_SCORE\]', response_text)
                if score_match:
                    text_analysis = response_text[:score_match.start()].strip()
                else:
                    text_analysis = response_text.strip()
            
            # extract vq_score - improved regex to handle various formats
            score_match = re.search(r'\[VQ_SCORE\]\s*(\d+)\s*\[/VQ_SCORE\]', response_text)
            if score_match:
                vq_score = int(score_match.group(1))
                # ensure score is within valid range
                vq_score = max(0, min(100, vq_score))
            else:
                # alternative pattern matching for score at the end
                score_match = re.search(r'(?:vq[_\s]*score|score)[:\s]*(\d+)', response_text, re.IGNORECASE)
                if score_match:
                    vq_score = int(score_match.group(1))
                    vq_score = max(0, min(100, vq_score))
                else:
                    # Additional fallback: look for any number in VQ_SCORE tags
                    score_match = re.search(r'\[VQ_SCORE\](.*?)\[/VQ_SCORE\]', response_text, re.DOTALL)
                    if score_match:
                        score_text = score_match.group(1).strip()
                        # Extract first number found
                        number_match = re.search(r'\d+', score_text)
                        if number_match:
                            vq_score = int(number_match.group())
                            vq_score = max(0, min(100, vq_score))
        
        except Exception as e:
            print(f"Error parsing LLM response: {e}")
            print(f"Response text: {response_text[:200]}...")
            # if parsing fails, return the full response as text analysis
            text_analysis = response_text.strip()
        
        return text_analysis, vq_score

    # process each viral hit
    if "Viral_Hits" in data:
        for i, viral_hit in enumerate(data["Viral_Hits"]):
            query_id = viral_hit.get("QueryID")
            
            # ai_summary field 
            if "AI_summary" not in viral_hit:
                viral_hit["AI_summary"] = ""
            
            # taxonomic information
            taxonomy_info = None
            
            # find the taxonomy information in this order: genus, family, order, norank
            for tax_level in ["Genus", "Family", "Order", "NoRank"]:
                if viral_hit.get(tax_level) and not pd.isna(viral_hit.get(tax_level)) and viral_hit.get(tax_level) != "":
                    # norank specially - remove "unclassified " prefix if present
                    search_name = viral_hit.get(tax_level)
                    if tax_level == "NoRank" and search_name.startswith("unclassified "):
                        search_name = search_name.replace("unclassified ", "", 1)
                    
                    # search tax level and name
                    matching_rows = taxonomy_df[(taxonomy_df['type'] == tax_level if tax_level != "NoRank" else 
                                            (taxonomy_df['type'].isin(['Genus', 'Family', 'Order']))) & 
                                            (taxonomy_df['name'] == search_name)]
                    
                    if not matching_rows.empty:
                        taxonomy_info = matching_rows.iloc[0]['info']
                        break
            
            # if found taxonomy info, generate ai summary and new vq_score
            if taxonomy_info:
                # prepare sequence data (exclude fullseq to avoid token limit issues)
                sequence_data = {k: v for k, v in viral_hit.items() if k != "FullSeq"}
                
                # related hmm hits for this queryid (exclude fullseq to avoid token limit issues)
                hmm_hits = []
                if "HMM_hits" in data:
                    hmm_hits = [{k: v for k, v in hit.items() if k != "FullSequence"} for hit in data["HMM_hits"] if hit.get("QueryID") == query_id]                
                
                # complete data for this sequence
                complete_data = {
                    "Viral_Hit": sequence_data,
                    "HMM_Hits": hmm_hits
                }
                
                # convert to string without the sequence data
                complete_data_str = json.dumps(complete_data, indent=2)
                
                prompt_input_data = {
                    "sequence_data": complete_data_str,
                    "taxonomy_info": taxonomy_info
                }
                
                # ai summary and vq_score based on model type
                try:
                    start_time = time.time()
                    
                    if args.model_type == "ollama":
                        user_content = ollama_user_prompt_template.format(**prompt_input_data)
                        
                        # ollama_options = {'temperature': 0.3}
                        
                        response_data = ollama.chat(
                            model=args.model_name,
                            messages=[
                                {'role': 'system', 'content': ollama_system_prompt},
                                {'role': 'user', 'content': user_content}
                            ],
                           
                        )
                        result = response_data['message']['content']  
                    
                    else:  # langchain models
                        if not langchain_llm_model:
                            print(f"CRITICAL: Langchain model for {args.model_type} not initialized for {query_id}.")
                            viral_hit["AI_summary"] = f"Error: Langchain model {args.model_type} not initialized."
                            continue  # skip to next viral hit
                        
                        # apply rate limiting
                        apply_api_rate_limit()
                    
                        # for chat models, use the chain with dictionary input
                        prompt = ChatPromptTemplate.from_template(langchain_template)
                        chain = prompt | langchain_llm_model
                        result = chain.invoke(prompt_input_data)
                    
                    end_time = time.time()
                    processing_time = end_time - start_time
                    
                    # remove thinking tags
                    cleaned_result = remove_thinking_tags(result)
                    
                    # parse the response
                    text_analysis, new_vq_score = parse_llm_response(cleaned_result)
                    
                    viral_hit["AI_summary"] = text_analysis.strip()
                    
                    # update vq_score
                    if new_vq_score is not None:
                        viral_hit["vq_score"] = new_vq_score
                    
                    
                except Exception as e:
                    print(f"Error generating summary for {query_id} with model {args.model_name}: {e}")
                    viral_hit["AI_summary"] = f"Error generating summary: {str(e)}"

            else:
                

                # prepare sequence data (exclude fullseq to avoid token limit issues)
                sequence_data = {k: v for k, v in viral_hit.items() if k != "FullSeq"}
                
                # related hmm hits for this queryid (exclude fullseq to avoid token limit issues)
                hmm_hits = []
                if "HMM_hits" in data:
                    hmm_hits = [{k: v for k, v in hit.items() if k != "FullSequence"} for hit in data["HMM_hits"] if hit.get("QueryID") == query_id]                
                
                # complete data for this sequence
                complete_data = {
                    "Viral_Hit": sequence_data,
                    "HMM_Hits": hmm_hits
                }
                
                # convert to string without the sequence data
                complete_data_str = json.dumps(complete_data, indent=2)
                
                prompt_input_data = {
                    "sequence_data": complete_data_str
                }
                
                # generate ai summary and vq_score based on model type
                try:
                    start_time = time.time()
                    
                    if args.model_type == "ollama":
                        user_content = ollama_user_prompt_template_noTax.format(**prompt_input_data)
                        
                        # ollama options (e.g., temperature)
                        # ollama_options = {'temperature': 0.3}
                        
                        response_data = ollama.chat(
                            model=args.model_name,
                            messages=[
                                {'role': 'system', 'content': ollama_system_prompt_noTax},
                                {'role': 'user', 'content': user_content}
                            ],
                            # options=ollama_options  # uncomment to use options
                        )
                        result = response_data['message']['content']  # string response
                    
                    else:  # langchain models (openai, anthropic, google)
                        if not langchain_llm_model:
                            print(f"CRITICAL: Langchain model for {args.model_type} not initialized for {query_id}.")
                            viral_hit["AI_summary"] = f"Error: Langchain model {args.model_type} not initialized."
                            continue  # skip to next viral hit
                        
                        # apply rate limiting for api calls only
                        apply_api_rate_limit()
                        
                        # for chat models, use the chain with dictionary input
                        prompt = ChatPromptTemplate.from_template(langchain_template_noTax)
                        chain = prompt | langchain_llm_model
                        result = chain.invoke(prompt_input_data)
                    
                    end_time = time.time()
                    processing_time = end_time - start_time
                    
                    # remove thinking tags if present
                    cleaned_result = remove_thinking_tags(result)
                    
                    # parse the response to extract text analysis and vq_score
                    text_analysis, new_vq_score = parse_llm_response(cleaned_result)
                    
                    viral_hit["AI_summary"] = text_analysis.strip()
                    
                    # update vq_score if a new score was provided
                    if new_vq_score is not None:
                        viral_hit["vq_score"] = new_vq_score
                        #print(f"Updated vq_score for {query_id}: {new_vq_score}")
                    
                    # print(f"Generated AI summary and vq_score for {query_id} in {processing_time:.2f} seconds.")
                    
                except Exception as e:
                    print(f"Error generating summary for {query_id} with model {args.model_name}: {e}")
                    viral_hit["AI_summary"] = f"Error generating summary: {str(e)}"


    # sort by biggest vq_score
    # first, normalize all vq_score values to int
    if 'Viral_Hits' in data:
        for hit in data['Viral_Hits']:
            score = hit.get('vq_score')
            if score is not None:
                hit['vq_score'] = int(score)  # convert string to int if needed
    
    # Now sort safely
    data['Viral_Hits'] = sorted(data['Viral_Hits'], key=lambda x: x['vq_score'], reverse=True)

    # save the updated json
    with open(output_json_path, 'w') as file:
        json.dump(data, file, indent=4)
    
    return data