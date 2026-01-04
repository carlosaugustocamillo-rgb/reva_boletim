import os
import json
import requests
import uuid
from datetime import datetime
import google.generativeai as genai
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure Gemini
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

# WhatsApp API Configuration (User will provide these later)
WA_PHONE_ID = os.environ.get("WA_PHONE_NUMBER_ID")
WA_TOKEN = os.environ.get("WA_ACCESS_TOKEN")
WA_BUSINESS_ID = os.environ.get("WA_BUSINESS_ACCOUNT_ID")
WA_API_VERSION = "v17.0"

def get_whatsapp_api_url():
    return f"https://graph.facebook.com/{WA_API_VERSION}/{WA_PHONE_ID}/messages"

def generate_draft_text(source_type, content_data):
    """
    Uses Gemini to generate a short, engaging WhatsApp message based on the content.
    
    Args:
        source_type (str): 'revacast_weekly' or 'revamais'
        content_data (dict): Contains 'title', 'summary', 'link', etc.
    
    Returns:
        str: The generated text for the WhatsApp message.
    """
    model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
    
    if source_type == 'revacast_weekly':
        prompt = f"""
        You are the social media manager for 'RevaCast', a scientific podcast for physiotherapists.
        Write a SHORT, EXCITING WhatsApp notification message (max 500 chars) announcing a new weekly episode.
        
        Content Highlights:
        {content_data.get('summary', 'New episode available!')}
        
        Link: {content_data.get('link', '[Link]')}
        
        Tone: Professional but energetic. Use emojis.
        Format:
        🔔 *RevaCast Weekly #[Date]*
        
        [Brief teaser of 1-2 key studies]
        
        🎧 Listen now: [Link]
        """
    elif source_type == 'revamais':
        prompt = f"""
        You are the editor of 'Reva+', a health newsletter for patients.
        Write a WARM, INVITING WhatsApp message (max 400 chars) announcing a new edition.
        
        Topic: {content_data.get('title')}
        Key Takeaway: {content_data.get('summary', 'Learn more about your health.')}
        Link: {content_data.get('link', '[Link]')}
        
        Tone: Empathetic, clear, educational. Use emojis.
        Format:
        🌿 *Reva+ News*
        
        [One sentence question to hook the reader about the topic]
        [Brief solution/insight]
        
        📖 Read more: [Link]
        """
    else:
        return "New content available!"

    try:
        response = model.generate_content(prompt)
        return response.text.strip()
    except Exception as e:
        print(f"⚠️ Error generating WhatsApp draft: {e}")
        return f"🔔 New {source_type} content available: {content_data.get('title')} - {content_data.get('link')}"

def create_draft(source_type, content_data):
    """
    Generates a draft and saves it to Firestore.
    """
    from firebase_service import save_firestore_document
    
    draft_text = generate_draft_text(source_type, content_data)
    draft_id = str(uuid.uuid4())
    
    draft_data = {
        "id": draft_id,
        "source": source_type,
        "content_ref": content_data.get('link'),
        "status": "DRAFT", # DRAFT, SCHEDULED, SENT
        "generated_text": draft_text,
        "original_data": content_data,
        "created_at": datetime.now().isoformat(),
        "scheduled_for": None
    }
    
    # Save to 'whatsapp_drafts' collection
    success = save_firestore_document("whatsapp_drafts", draft_id, draft_data)
    
    if success:
        print(f"✅ WhatsApp Draft created: {draft_id}")
        return draft_data
    else:
        print("❌ Failed to save WhatsApp Draft.")
        return None

def send_message(to_number, text_body):
    """
    Sends a text message using WhatsApp Cloud API.
    """
    if not (WA_PHONE_ID and WA_TOKEN):
        return {"error": "Missing WhatsApp Credentials (WA_PHONE_NUMBER_ID or WA_ACCESS_TOKEN)"}

    url = get_whatsapp_api_url()
    headers = {
        "Authorization": f"Bearer {WA_TOKEN}",
        "Content-Type": "application/json"
    }
    
    payload = {
        "messaging_product": "whatsapp",
        "to": to_number,
        "type": "text",
        "text": {"body": text_body}
    }
    
    try:
        r = requests.post(url, headers=headers, json=payload)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return {"error": str(e), "details": r.text if 'r' in locals() else "Request failed"}

def send_template_message(to_number, template_name="hello_world", lang_code="en_US"):
    """
    Sends a template message (Required for starting conversations).
    """
    if not (WA_PHONE_ID and WA_TOKEN):
        return {"error": "Missing WhatsApp Credentials"}

    url = get_whatsapp_api_url()
    headers = {
        "Authorization": f"Bearer {WA_TOKEN}",
        "Content-Type": "application/json"
    }
    
    payload = {
        "messaging_product": "whatsapp",
        "to": to_number,
        "type": "template",
        "template": {
            "name": template_name,
            "language": {
                "code": lang_code
            }
        }
    }
    
    try:
        r = requests.post(url, headers=headers, json=payload)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return {"error": str(e), "details": r.text if 'r' in locals() else "Request failed"}


def list_drafts(status_filter=None):
    """
    Lists drafts from Firestore.
    """
    from firebase_service import get_firestore_db
    
    db = get_firestore_db()
    if not db:
        print("❌ Firestore not available for listing drafts.")
        return []
    
    try:
        collection_ref = db.collection("whatsapp_drafts")
        if status_filter:
            query = collection_ref.where("status", "==", status_filter)
        else:
            # Order by created_at desc by default if possible, but requires index.
            # For now just get all.
            query = collection_ref
            
        docs = query.stream()
        return [doc.to_dict() for doc in docs]
    except Exception as e:
        print(f"Error listing drafts: {e}")
        return []

def update_draft_status(draft_id, status, scheduled_date=None):
    """
    Updates the status of a draft.
    """
    from firebase_service import update_firestore_document
    
    update_data = {"status": status}
    if scheduled_date:
        update_data["scheduled_for"] = scheduled_date
        
    return update_firestore_document("whatsapp_drafts", draft_id, update_data)

