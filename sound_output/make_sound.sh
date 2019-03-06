# Copyright (c) 2019, Alexander D. Kaiser
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


ffmpeg -i mitral_1202357_512_git_3735567_two_leaflet_eight_connectors_real_time.mp4 \
    -i heart_sounds_and_silence.wav -c:v libx264 -strict -2 -map 0:v -map 1:a mitral_1202357_512_git_3735567_two_leaflet_eight_connectors_real_time_SOUND.mp4


#ffmpeg -i mitral_1092466_256_git_4dc23b8_two_leaflet_updated_damping_damping_2x_real_time.mp4 \
#    -i heart_sounds_and_silence.wav -c:v libx264 -strict -2 -map 0:v -map 1:a mitral_1092466_256_git_4dc23b8_two_leaflet_updated_damping_damping_2x_real_time_SOUND.mp4


# -c:a libvorbis
# -af "volume=10.0"

# works great, some annoying clicks
# ffmpeg -i mitral_cycle_PERIODIC_9004892_git_d7914ba6_more_wrap_256_radial_real_time.mp4 \
#        -t 2.3 -i heart_sounds.wav -c:v libx264 -strict -2 -map 0:v -map 1:a mitral_cycle_PERIODIC_9004892_git_d7914ba6_more_wrap_256_radial_real_time_SOUND.mp4

