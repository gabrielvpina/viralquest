import json
import pandas as pd
import os
import re
import argparse
import time
import ollama # Adicionado para uso direto do Ollama

# Langchain imports são condicionais ou para tipos de modelo específicos
from langchain_core.prompts import ChatPromptTemplate

def remove_thinking_tags(text):
    """
    Remove <thinking>...</thinking> tags and their content from the text.
    """
    if hasattr(text, 'content'):
        text_content = text.content
    elif isinstance(text, dict) and 'content' in text:
        text_content = text['content']
    else:
        text_content = str(text)
    
    cleaned_text = re.sub(r'<thinking>.*?</thinking>', '', text_content, flags=re.DOTALL)
    cleaned_text = re.sub(r'<think>.*?</think>', '', cleaned_text, flags=re.DOTALL)
    cleaned_text = re.sub(r'<[Tt]hinking>.*?</[Tt]hinking>', '', cleaned_text, flags=re.DOTALL)
    cleaned_text = re.sub(r'<[Tt]hink>.*?</[Tt]hink>', '', cleaned_text, flags=re.DOTALL)
    cleaned_text = re.sub(r'\s+', ' ', cleaned_text).strip()
    return cleaned_text

def get_llm_model(args):
    """
    Initialize and return the appropriate LLM based on user arguments.
    For 'ollama' type, this function will return None as it's handled directly.
    """
    if args.model_type == "ollama":
        # Não é necessário um objeto LLM Langchain para uso direto do ollama
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
    json_file_path = args.input_json
    taxonomy_data_path = args.taxonomy_data
    output_json_path = args.output_json if args.output_json else json_file_path
    
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    taxonomy_df = pd.read_csv(taxonomy_data_path, sep='\t')
    
    # Inicializar o modelo LLM Langchain (apenas se não for Ollama local)
    langchain_llm_model = None
    if args.model_type != "ollama":
        try:
            langchain_llm_model = get_llm_model(args)
            if langchain_llm_model is None: # Deveria ter sido tratado por get_llm_model, mas checagem dupla
                 raise ValueError(f"Model could not be initialized for type: {args.model_type}")
            print(f"Successfully loaded {args.model_type} model: {args.model_name}")
        except Exception as e:
            print(f"Error initializing Langchain model: {e}")
            return None
    elif args.model_type == "ollama":
        try:
            # Verifica se o servidor Ollama está acessível e o modelo existe
            ollama.show(args.model_name) 
            print(f"Successfully connected to Ollama. Using model: {args.model_name}")
        except Exception as e:
            print(f"Error connecting to Ollama or model '{args.model_name}' not found/available: {e}")
            print(f"Please ensure Ollama is running and the model is downloaded (e.g., 'ollama pull {args.model_name}').")
            return None

    # Template para modelos Langchain (OpenAI, Anthropic, Google)
    langchain_template_string = """
    You are an expert in virology and will analyze some information of potential viral sequences from a bioinformatics analysis.
    You need to be skeptical in your analysis. 
    Here is the information:
    
    {sequence_data}
    
    Here is the information about the viral order/family/genus:
    {taxonomy_info}
    
    Based on this information, make a text (maximum of 200 words) about the taxonomy information (hosts), 
    and evaluate if the sequence above is a known virus, a novel virus, or a non-viral sequence. 
    Talk about the demarcation criteria of the taxonomy group. 
    Remember, the analyses also depends of the results in BLASTx and BLASTn identity and Coverage, 
    normally a known virus has more than 90% in aminoacid/nucleotide similarity and at least 70% of coverage.
    """

    # Prompts para uso direto do Ollama
    ollama_system_prompt = """You are an expert in virology and will analyze some information of potential viral sequences from a bioinformatics analysis.
You need to be skeptical in your analysis. 
Based on the information provided by the user, make a text (maximum of 200 words) about the taxonomy information (hosts), 
and evaluate if the sequence is a known virus, a novel virus, or a non-viral sequence. 
Talk about the demarcation criteria of the taxonomy group. 
Remember, the analyses also depends of the results in BLASTx and BLASTn identity and Coverage, 
normally a known virus has more than 90% in aminoacid/nucleotide similarity and at least 70% of coverage."""

    ollama_user_prompt_template = """Here is the sequence information:
{sequence_data}

Here is the information about the viral order/family/genus:
{taxonomy_info}"""
    
    if "Viral_Hits" in data:
        for i, viral_hit in enumerate(data["Viral_Hits"]):
            query_id = viral_hit.get("QueryID")
            
            if "AI_summary" not in viral_hit:
                viral_hit["AI_summary"] = ""
            
            taxonomy_info = None
            for tax_level in ["Genus", "Family", "Order", "NoRank"]:
                if viral_hit.get(tax_level) and not pd.isna(viral_hit.get(tax_level)) and viral_hit.get(tax_level) != "":
                    search_name = viral_hit.get(tax_level)
                    if tax_level == "NoRank" and search_name.startswith("unclassified "):
                        search_name = search_name.replace("unclassified ", "", 1)
                    
                    matching_rows = taxonomy_df[(taxonomy_df['type'] == tax_level if tax_level != "NoRank" else 
                                              (taxonomy_df['type'].isin(['Genus', 'Family', 'Order']))) & 
                                              (taxonomy_df['name'] == search_name)]
                    
                    if not matching_rows.empty:
                        taxonomy_info = matching_rows.iloc[0]['info']
                        print(f"Found taxonomy info for {query_id} using {tax_level}: {search_name}")
                        break
            
            if taxonomy_info:
                sequence_data = {k: v for k, v in viral_hit.items() if k != "FullSeq"}
                hmm_hits = []
                if "HMM_hits" in data:
                    hmm_hits = [{k: v for k, v in hit.items() if k != "FullSequence"} for hit in data["HMM_hits"] if hit.get("QueryID") == query_id]
                
                complete_data = {
                    "Viral_Hit": sequence_data,
                    "HMM_Hits": hmm_hits
                }
                complete_data_str = json.dumps(complete_data, indent=2)
                
                prompt_input_data = {
                    "sequence_data": complete_data_str,
                    "taxonomy_info": taxonomy_info
                }

                raw_llm_response = None
                try:
                    start_time = time.time()
                    
                    if args.model_type == "ollama":
                        user_content = ollama_user_prompt_template.format(**prompt_input_data)
                        
                        print(f"\n--- System Prompt for {query_id} (Ollama Direct) ---")
                        print(ollama_system_prompt)
                        print(f"--- User Prompt for {query_id} (Ollama Direct) ---")
                        print(user_content)
                        print("--- End of Prompts ---")

                        # Opções para o Ollama (ex: temperatura)
                        # ollama_options = {'temperature': 0.3}

                        response_data = ollama.chat(
                            model=args.model_name,
                            messages=[
                                {'role': 'system', 'content': ollama_system_prompt},
                                {'role': 'user', 'content': user_content}
                            ],
                            # options=ollama_options # Descomente para usar opções
                        )
                        raw_llm_response = response_data['message']['content'] # String da resposta
                    
                    else: # Modelos Langchain (OpenAI, Anthropic, Google)
                        if not langchain_llm_model:
                            print(f"CRITICAL: Langchain model for {args.model_type} not initialized for {query_id}.")
                            viral_hit["AI_summary"] = f"Error: Langchain model {args.model_type} not initialized."
                            continue # Pula para o próximo hit viral

                        prompt_template_obj = ChatPromptTemplate.from_template(langchain_template_string)
                        prompt_messages = prompt_template_obj.format_messages(**prompt_input_data)
                        full_prompt_text = "\n".join([msg.content for msg in prompt_messages if hasattr(msg, 'content')])
                        
                        print(f"\n--- Prompt for {query_id} ({args.model_type} via Langchain) ---")
                        print(full_prompt_text)
                        print("--- End of Prompt ---")
                        
                        chain = prompt_template_obj | langchain_llm_model
                        raw_llm_response = chain.invoke(prompt_input_data) # Objeto de mensagem Langchain
                    
                    end_time = time.time()
                    processing_time = end_time - start_time
                    
                    cleaned_result = remove_thinking_tags(raw_llm_response)
                    viral_hit["AI_summary"] = cleaned_result.strip()
                    print(f"Generated AI summary for {query_id} in {processing_time:.2f} seconds.")
                
                except Exception as e:
                    print(f"Error generating summary for {query_id} with model {args.model_name}: {e}")
                    viral_hit["AI_summary"] = f"Error generating summary: {str(e)}"
            else:
                print(f"No matching taxonomy information found for {query_id}")
    
    with open(output_json_path, 'w') as file:
        json.dump(data, file, indent=4)
    
    print(f"Updated JSON saved to {output_json_path}")
    return data

def main():
    parser = argparse.ArgumentParser(description='Analyze viral sequences using AI models')
    parser.add_argument('--input-json', required=True, help='Path to the input JSON file with viral sequence data')
    parser.add_argument('--taxonomy-data', required=True, help='Path to the CSV file with viral taxonomy information')
    parser.add_argument('--model-type', required=True, choices=['ollama', 'openai', 'anthropic', 'google'], 
                        help='Type of model to use for analysis (ollama for local direct usage)')
    parser.add_argument('--model-name', required=True, 
                        help='Name of the model (e.g., "qwen2:1.5b" for ollama, "gpt-3.5-turbo" for OpenAI)')
    parser.add_argument('--api-key', help='API key for cloud models (OpenAI, Anthropic, Google)')
    parser.add_argument('--output-json', help='Path to save the output JSON (defaults to overwriting input)')
    args = parser.parse_args()
    analyze_viral_sequences(args)

if __name__ == "__main__":
    main()