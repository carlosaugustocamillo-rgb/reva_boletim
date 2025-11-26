import google.generativeai as genai
import os
from dotenv import load_dotenv

load_dotenv()
genai.configure(api_key=os.getenv("GEMINI_API_KEY"))

def test_image_gen():
    model_name = "imagen-4.0-generate-001"
    print(f"Testing {model_name}...")
    try:
        model = genai.GenerativeModel(model_name)
        response = model.generate_content("Generate an image of a futuristic hospital.")
        
        print("Response parts:", len(response.parts))
        for part in response.parts:
            print("Part:", part)
            if hasattr(part, 'inline_data') and part.inline_data:
                print("Found inline data!")
                # Save it
                import base64
                img_data = part.inline_data.data
                # img_data is bytes? or base64 string? usually bytes in the object
                # actually part.inline_data.data is bytes
                with open("test_gemini_image.png", "wb") as f:
                    f.write(img_data)
                print("Saved test_gemini_image.png")
                return
            
            # Check for other formats
            print(dir(part))
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_image_gen()
