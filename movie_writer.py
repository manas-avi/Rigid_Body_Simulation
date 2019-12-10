import cv2
import os

image_folder = 'imgs'
video_name = 'video.avi'

# images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
img = os.listdir(image_folder)[0]
frame = cv2.imread(os.path.join(image_folder, img))
height, width, layers = frame.shape

video = cv2.VideoWriter(video_name, 0, 30, (width,height))

num_files = len(os.listdir(image_folder))
for count in range(1,num_files+1):
	img = cv2.imread(os.path.join(image_folder, str(count) + '.png'))
	video.write(img)

cv2.destroyAllWindows()
video.release()
