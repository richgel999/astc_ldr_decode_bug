#define TEST_LDR_SRGB_DECODE (0)

// basisu_tool.cpp
// Copyright (C) 2019-2022 Binomial LLC. All Rights Reserved.
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#if _MSC_VER
// For sprintf(), strcpy() 
#define _CRT_SECURE_NO_WARNINGS (1)
#endif

#define ASTC_HELPERS_IMPLEMENTATION (1)
#include "astc_helpers.h"

#include "transcoder/basisu.h"
#include "transcoder/basisu_transcoder_internal.h"
#include "encoder/basisu_enc.h"
#include "encoder/basisu_etc.h"
#include "encoder/basisu_gpu_texture.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

using namespace basisu;

#include "astc-encoder/astcenc.h"
#include "encoder/basisu_astc_decomp.h"

struct astc_blk
{
	uint8_t m_vals[16];
};

astcenc_context* g_astc_context;

static void decomp_astc_block_init(uint32_t w, uint32_t h, astcenc_profile profile = ASTCENC_PRF_HDR)
{
	astcenc_config config;

	astcenc_error status = astcenc_config_init(profile, w, h, 1, ASTCENC_PRE_MEDIUM, ASTCENC_FLG_DECOMPRESS_ONLY, &config);
	if (status != ASTCENC_SUCCESS)
	{
		printf("astcenc_config_init() failed\n");
		exit(1);
	}

	status = astcenc_context_alloc(&config, 1, &g_astc_context);
	if (status != ASTCENC_SUCCESS)
	{
		printf("astcenc_context_alloc() failed\n");
		exit(1);
	}
}

static void decomp_astc_block_deinit()
{
	astcenc_context_free(g_astc_context);
}

static bool decomp_astc_block_f32(const astc_blk& blk, vec4F* pPixels, uint32_t block_w, uint32_t block_h)
{
	static const astcenc_swizzle swizzle{
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_image astc_image;
	astc_image.dim_x = block_w;
	astc_image.dim_y = block_h;
	astc_image.dim_z = 1;
	astc_image.data_type = ASTCENC_TYPE_F32;

	uint8_t* slices = (uint8_t*)pPixels;
	astc_image.data = reinterpret_cast<void**>(&slices);

	int status = astcenc_decompress_image(g_astc_context, (const uint8_t*)&blk, 16, &astc_image, &swizzle, 0);

	return (status == ASTCENC_SUCCESS);
}

static bool decomp_astc_block_u8(const astc_blk& blk, uint8_t* pPixels, uint32_t block_w, uint32_t block_h)
{
	static const astcenc_swizzle swizzle{
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_image astc_image;
	astc_image.dim_x = block_w;
	astc_image.dim_y = block_h;
	astc_image.dim_z = 1;
	astc_image.data_type = ASTCENC_TYPE_U8;

	uint8_t* slices = (uint8_t*)pPixels;
	astc_image.data = reinterpret_cast<void**>(&slices);

	int status = astcenc_decompress_image(g_astc_context, (const uint8_t*)&blk, 16, &astc_image, &swizzle, 0);

	return (status == ASTCENC_SUCCESS);
}

static int test()
{
	astc_helpers::init_tables();

	const uint32_t block_w = 12, block_h = 12;
	
	if (!astc_helpers::is_valid_block_size(block_w, block_h))
	{
		printf("Bad block size\n");
		exit(1);
	}

	int min_endpoint_ise_range = -1;
	const bool ldr_cems_only = true;
	const bool hdr_cems_only = false;
	const uint32_t first_cem = 0;
	const uint32_t last_cem = 15;

	bool test_srgb = ldr_cems_only && (TEST_LDR_SRGB_DECODE != 0);
	decomp_astc_block_init(block_w, block_h, ldr_cems_only ? (test_srgb ? ASTCENC_PRF_LDR_SRGB : ASTCENC_PRF_LDR) : ASTCENC_PRF_HDR);

	uint32_t total_modes = 0;

	basisu::rand rnd;
	rnd.seed(1000 + block_w * 52772 + block_h * 10055);

	for (uint32_t dp = 0; dp < 2; dp++)
	{
		for (uint32_t num_parts = 1; num_parts <= (dp ? 3 : 4); num_parts++)
		{
			for (uint32_t weight_range = 0; weight_range <= 11; weight_range++)
			{
				for (uint32_t cem0 = 0; cem0 <= 15; cem0++)
				{
					for (uint32_t cem1 = (num_parts >= 2) ? first_cem : cem0; cem1 <= ((num_parts >= 2) ? last_cem : cem0); cem1++)
					{
						for (uint32_t cem2 = (num_parts >= 3) ? first_cem : cem1; cem2 <= ((num_parts >= 3) ? last_cem : cem1); cem2++)
						{
							for (uint32_t cem3 = (num_parts >= 4) ? first_cem : cem2; cem3 <= ((num_parts == 4) ? last_cem : cem2); cem3++)
							{
								for (uint32_t w = 2; w <= block_w; w++)
								{
									for (uint32_t h = 2; h <= block_h; h++)
									{
										const bool is_ldr_cem0 = astc_helpers::is_cem_ldr(cem0);
										const bool is_ldr_cem1 = astc_helpers::is_cem_ldr(cem1);
										const bool is_ldr_cem2 = astc_helpers::is_cem_ldr(cem2);
										const bool is_ldr_cem3 = astc_helpers::is_cem_ldr(cem3);

										if (ldr_cems_only)
										{
											if ((!is_ldr_cem0) || (!is_ldr_cem1) || (!is_ldr_cem2) || (!is_ldr_cem3))
												continue;
										}

										if (hdr_cems_only)
										{
											if ((is_ldr_cem0) || (is_ldr_cem1) || (is_ldr_cem2) || (is_ldr_cem3))
												continue;
										}

										astc_helpers::log_astc_block log_block;
										log_block.clear();

										log_block.m_grid_width = w;
										log_block.m_grid_height = h;
										log_block.m_dual_plane = dp != 0;
										log_block.m_color_component_selector = 0;
										log_block.m_endpoint_ise_range = 20;
										log_block.m_weight_ise_range = weight_range;
										log_block.m_num_partitions = num_parts;
										log_block.m_partition_id = 0;
										log_block.m_color_endpoint_modes[0] = cem0;
										log_block.m_color_endpoint_modes[1] = (num_parts >= 2) ? cem1 : 0;
										log_block.m_color_endpoint_modes[2] = (num_parts >= 3) ? cem2 : 0;
										log_block.m_color_endpoint_modes[3] = (num_parts >= 4) ? cem3 : 0;

										uint32_t total_endpoint_vals = 0;
										for (uint32_t i = 0; i < num_parts; i++)
											total_endpoint_vals += astc_helpers::get_num_cem_values(log_block.m_color_endpoint_modes[i]);
										if (total_endpoint_vals > astc_helpers::MAX_ENDPOINTS)
											continue;

										// Random endpoints
										uint32_t endpoint_levels = astc_helpers::get_ise_levels(log_block.m_endpoint_ise_range);
										for (uint32_t j = 0; j < total_endpoint_vals; j++)
											log_block.m_endpoints[j] = (uint8_t)rnd.irand(0, endpoint_levels - 1);
										
										const uint32_t weight_levels = astc_helpers::get_ise_levels(log_block.m_weight_ise_range);
										const uint32_t total_weights = (log_block.m_dual_plane ? 2 : 1) * log_block.m_grid_width * log_block.m_grid_height;
										if (total_weights > astc_helpers::MAX_GRID_WEIGHTS)
											continue;

										// Random weights
										for (uint32_t j = 0; j < total_weights; j++)
											log_block.m_weights[j] = (uint8_t)rnd.irand(0, weight_levels - 1);

										// Pack the logical block to a physical ASTC block
										astc_helpers::astc_block phys_block;
										int expected_endpoint_range = -1;
										bool status = astc_helpers::pack_astc_block(phys_block, log_block, &expected_endpoint_range);

										if (!status)
										{
											// Packed failed at the highest endpoint range, so try at the highest supported endpoint range for this config
											if (expected_endpoint_range != -1)
											{
												log_block.m_endpoint_ise_range = expected_endpoint_range;

												endpoint_levels = astc_helpers::get_ise_levels(log_block.m_endpoint_ise_range);
												for (uint32_t j = 0; j < total_endpoint_vals; j++)
													log_block.m_endpoints[j] = (uint8_t)rnd.irand(0, endpoint_levels - 1);

												status = astc_helpers::pack_astc_block(phys_block, log_block, &expected_endpoint_range);
											}

											if (!status)
												continue;
										}

										if ((int)log_block.m_endpoint_ise_range < min_endpoint_ise_range)
											continue;

										astc_helpers::log_astc_block unpacked_log_blk;

										// Now unpack the physical block back to logical and make sure it round-tripped correctly.
										{
											bool q = astc_helpers::unpack_block(&phys_block, unpacked_log_blk, block_w, block_h);
											if (!q)
											{
												printf("Failure\n");
												exit(1);
											}

											bool same = log_block.m_grid_width == unpacked_log_blk.m_grid_width;
											same = same & (log_block.m_grid_height == unpacked_log_blk.m_grid_height);
											same = same & (log_block.m_dual_plane == unpacked_log_blk.m_dual_plane);
											same = same & (log_block.m_color_component_selector == unpacked_log_blk.m_color_component_selector);
											same = same & (log_block.m_endpoint_ise_range == unpacked_log_blk.m_endpoint_ise_range);
											same = same & (log_block.m_weight_ise_range == unpacked_log_blk.m_weight_ise_range);
											same = same & (log_block.m_num_partitions == unpacked_log_blk.m_num_partitions);
											same = same & (log_block.m_partition_id == unpacked_log_blk.m_partition_id);
											for (uint32_t i = 0; i < log_block.m_num_partitions; i++)
											{
												same = same & (log_block.m_color_endpoint_modes[i] == unpacked_log_blk.m_color_endpoint_modes[i]);
											}
											for (uint32_t i = 0; i < 64; i++)
											{
												same = same & (log_block.m_weights[i] == unpacked_log_blk.m_weights[i]);
											}
											for (uint32_t i = 0; i < 18; i++)
											{
												same = same & (log_block.m_endpoints[i] == unpacked_log_blk.m_endpoints[i]);
											}
											if (!same)
											{
												printf("Failure\n");
												exit(1);
											}

										}

										printf("%u: blk: %ux%u, cems: %u %u %u %u, grid: %ux%u, dp: %u, parts: %u, endpoint range: %i, endpoint levels: %i, weight range: %u, weight levels: %u\n",
											total_modes,
											block_w, block_h,
											cem0, cem1, cem2, cem3, w, h, dp, num_parts,
											log_block.m_endpoint_ise_range, astc_helpers::get_ise_levels(log_block.m_endpoint_ise_range),
											weight_range, astc_helpers::get_ise_levels(weight_range));
										total_modes++;

										printf("Physical ASTC block bytes: { ");
										for (uint32_t i = 0; i < 16; i++)
											printf("0x%X%c", ((uint8_t*)&phys_block)[i], (i != 15) ? ',' : ' ');
										printf("} \n");

										// Now decompress the block to pixels

										if (!is_ldr_cem0 || !is_ldr_cem1 || !is_ldr_cem2 || !is_ldr_cem3)
										{
											// The block has 1 or more HDR CEMs
											
											// Google's decompressor
											vec4F pixels[12 * 12];

											status = basisu_astc::astc::decompress_hdr(&pixels[0][0], (uint8_t*)&phys_block, block_w, block_h);
											if (!status)
											{
												assert(0);
												continue;
											}

											// My decompressor
											uint16_t alt_pixels[12 * 12 * 4];
											status = astc_helpers::decode_block(unpacked_log_blk, alt_pixels, block_w, block_h, astc_helpers::cDecodeModeHDR16);
											if (!status)
											{
												assert(0);
												continue;
											}
											
											// ARM's decompressor
											vec4F arm_pixels[12 * 12];
											status = decomp_astc_block_f32((astc_blk&)phys_block, arm_pixels, block_w, block_h);
											if (!status)
											{
												printf("decomp_astc_block() failed!\n");
												exit(1);
											}

											// Compare the outputs
											for (uint32_t y = 0; y < block_h; y++)
											{
												for (uint32_t x = 0; x < block_w; x++)
												{
													for (uint32_t c = 0; c < 4; c++)
													{
														float a = pixels[x + y * block_w][c];
														
														half_float hb = alt_pixels[(x + y * block_w) * 4 + c];
														float b = half_to_float(hb);

														float d = arm_pixels[x + y * block_w][c];

														if ((a != b) || (a != d))
														{
															printf("Pixel compare failure %ux%u comp %u\n", x, y, c);
														}
													}
												}
											}

										}
										else
										{
											// The block uses all LDR CEM's
											
											// Google's decompressor
											uint8_t pixels[12 * 12 * 4];
											status = basisu_astc::astc::decompress_ldr(pixels, (uint8_t*)&phys_block, test_srgb, block_w, block_h);
											if (!status)
											{
												assert(0);
												continue;
											}

											// My decompressor
											uint8_t alt_pixels[12 * 12 * 4];
											status = astc_helpers::decode_block(unpacked_log_blk, alt_pixels, block_w, block_h, test_srgb ? astc_helpers::cDecodeModeSRGB8 : astc_helpers::cDecodeModeLDR8);
											if (!status)
											{
												assert(0);
												continue;
											}

											// ARM's decompressor
											uint8_t arm_pixels[12 * 12 * 4];
											status = decomp_astc_block_u8((astc_blk&)phys_block, arm_pixels, block_w, block_h);
											if (!status)
											{
												printf("decomp_astc_block() failed!\n");
												exit(1);
											}

											// Compare the outputs
											for (uint32_t y = 0; y < block_h; y++)
											{
												for (uint32_t x = 0; x < block_w; x++)
												{
													for (uint32_t c = 0; c < 4; c++)
													{
														if ((pixels[(x + y * block_w) * 4 + c] != alt_pixels[(x + y * block_w) * 4 + c]) || 
															//(iabs(pixels[(x + y * block_w) * 4 + c] - arm_pixels[(x + y * block_w) * 4 + c]) > 1) 
															(pixels[(x + y * block_w) * 4 + c] != arm_pixels[(x + y * block_w) * 4 + c])
														    )
														{
															printf("Pixel compare failure %ux%u comp %u, Google: %u, mine: %u, ARM: %u\n", x, y, c,
																pixels[(x + y * block_w) * 4 + c],
																alt_pixels[(x + y * block_w) * 4 + c],
																arm_pixels[(x + y * block_w) * 4 + c]);
														}
													}
												}
											}
										}
																				

									} // h

								} // w

							} // cem 3
						} // cem 2
					} // cem1				
				} // cem0

			} // weight_range

		} // num_parts
	} // dp

	printf("OK\n");

	return EXIT_SUCCESS;
}

int main(int argc, const char** argv)
{
	basisu_encoder_init(false, false);

	return test();
}

