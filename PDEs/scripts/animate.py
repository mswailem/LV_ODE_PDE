import os
from moviepy.video.io.ImageSequenceClip import ImageSequenceClip
import re  # For regular expression processing

# Custom sorting key function to extract numeric parts
def extract_number(filename):
    match = re.search(r"(\d+)", filename)  # Find the first numeric substring
    return int(match.group(1)) if match else 0  # Convert it to an integer, or return 0 if no match

# Get all PNG image files and sort them numerically
image_files = [file for file in sorted(os.listdir('.'), key=extract_number) if file.endswith('.png')]

# Create a video clip from the image sequence
video_clip = ImageSequenceClip(image_files, fps=10)  # Adjust fps to control speed

# Save the video
video_file = "animation.mp4"
video_clip.write_videofile(video_file, codec="libx264")

print(f"MP4 video saved to {video_file}")

