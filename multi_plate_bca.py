metadata = {
    "protocolName": "Multi-plate BCA protocol",
    "author": "Nico To (modification of Sasha's original BCA protocol)",
    "description": "BCA for 25+ samples",
}
requirements = {"robotType": "Flex", "apiLevel": "2.20"}
import math
from opentrons import protocol_api

def get_vol_50ml_falcon(height):
    """
    Get's the volume of the liquid in the tube
    Height: height of liquid in tube in mm (start from tube bottom)
    Return: volume of liquid in tube in µl
    """
    volume = (1000 * (height - 9)) / 1.8
    return volume
def get_height_50ml_falcon(volume):
    """
    Get's the height of the liquid in the tube
    Volume: volume of liquid in tube in µl
    Return: hieght from bottom of tube in millimeters
    """
    height = (1.8 * (volume / 1000)) + 9
    return height
def get_height_15ml_falcon(volume):
    """
    Get's the height of the liquid in the tube
    Volume: volume of liquid in tube in ul
    Return: height in mm from the bottom of tube that pipette should go to
    """
    volume = volume / 1000
    if volume <= 1:  # cone part
        # print(-3.33*(volume**2)+15.45*volume+9.50)
        return -3.33 * (volume**2) + 15.45 * volume + 9.50 - 1  # −3.33x2+15.45x+9.50
    else:
        return 6.41667 * volume + 15.1667 - 5
def get_vol_15ml_falcon(height):
    """
    Get's the volume of the liquid in the tube
    Height: height in mm from the bottom of tube that pipette should go to
    Return: volume of liquid in tube in ul
    """
    if height <= 20.62:  # cone part
        volume = (((15.45 + math.sqrt(351.9225 - (13.32 * height)))) / 6.66) * 1000
        return volume
    else:
        volume = ((height - 10.1667) / 6.41667) * 1000
        return volume

def add_parameters(parameters):
    parameters.add_int(
        variable_name="number_samples",
        display_name="number_samples",
        description="Number of input samples.",
        default=96,
        minimum=1,  #change to 25
        maximum=96,
        unit="samples",
    )
    parameters.add_bool(
        variable_name="dulute_with_walt",
        display_name="Dilute using Walter",
        description="True: Dilute using Walter, False: Dilute manually",
        default=True,
    )
    parameters.add_int(
        variable_name="sample_vol",
        display_name="Amount of Sample",
        description="Amount of sample required in the duilution (only required if dilute_with_walt is True)",
        default=10,
        minimum=5,
        maximum=200,
        unit="ul",
    )
    parameters.add_int(
        variable_name="buffer_vol",
        display_name="Amount of Buffer",
        description="Amount of buffer required in the duilution (only required if dilute_with_walt is True)",
        default=90,
        minimum=5,
        maximum=200,
        unit="ul",
    )
    parameters.add_int(
        variable_name="working_sample_vol",
        display_name="Working Sample Volume",
        description="Volume of working sample (Volume for the samples in the flat well plate)",
        default=25,
        minimum=5,
        maximum=50,
        unit="ul",
    )
    parameters.add_int(
        variable_name="tip_type",
        display_name="types of tips",
        choices=[
            {"display_name": "1000ul", "value": 1000},
            {"display_name": "200ul", "value": 200},
        ],
        default=200,
    )
    parameters.add_int(
        variable_name="replication_mode",
        display_name="Replication Mode",
        choices=[
            {"display_name": "triplicate", "value": 3},
            {"display_name": "duplicate", "value": 2},
        ],
        default=3,
    )

    parameters.add_bool(
        variable_name="dry_run",
        display_name="Dry Run",
        description="Return tips (ignore this unless you are testing)",
        default=False,
    )


def run(protocol: protocol_api.ProtocolContext):
    replication_mode= protocol.params.replication_mode
    number_samples = protocol.params.number_samples
    is_dry_run = protocol.params.dry_run
    working_sample_vol = protocol.params.working_sample_vol
    pipette_max = 200-5

    # LOADING TIPS
    tips = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_1000uL", slot)
        for slot in ["A3"]
    ]
    chute = protocol.load_waste_chute()

    def remove_tip(pipette):
        if is_dry_run:
            pipette.return_tip()
        else:
            pipette.drop_tip(chute)

    def pick_up(pip):
        nonlocal tips
        nonlocal staging_racks
        nonlocal count

        try:
            pip.tip_racks = tips
            pip.pick_up_tip()

        except protocol_api.labware.OutOfTipsError:
            check_tips()
            pick_up(pip)    
    def check_tips():
        nonlocal tips
        nonlocal staging_racks
        nonlocal count
        # tip_box = protocol.load_labware('opentrons_flex_96_filtertiprack_1000uL', 'A3')
        for i in range (0,1):
            tip_box_slots = ['A3']
            bottom_right_well = tips[i].wells_by_name()['H12']
            
            if bottom_right_well.has_tip or protocol.deck['D4'] == None:
                protocol.comment("A tip is present in the bottom-right corner (H12). or all staging slots are empty")
                if protocol.deck['D4'] == None:
                    protocol.comment("No tip box detected in slot D4.")
                    staging_slots = ['A4', 'B4', 'C4', 'D4']
                    staging_racks = [protocol.load_labware('opentrons_flex_96_filtertiprack_1000uL',
                                      slot) for slot in staging_slots]
                pass
            else:
                protocol.comment("\n\n\n Starting moving phase")
                protocol.move_labware(
                        labware=tips[i],
                        new_location=chute,
                        use_gripper=True
                    )
                rack_num = 0
                for slot in ['A4', 'B4', 'C4', 'D4']:
                    labware = protocol.deck[slot]
                    if labware and labware.is_tiprack:
                        tips[i] = staging_racks[rack_num]
                        protocol.move_labware(
                            labware=staging_racks[rack_num],
                            new_location=tip_box_slots[i],
                            use_gripper=True
                        )
                        break
                        # protocol.comment(f"A tip box is present in slot {slot}.")
                    else:
                        protocol.comment(f"No tip box detected in slot {slot}.")
                        rack_num+=1
                        pass

    def find_aspirate_height(pip, source_well):
        lld_height = (
            pip.measure_liquid_height(source_well) - source_well.bottom().point.z
        )
        aspirate_height = max(lld_height - 5, 1)
        return aspirate_height

    # LOADING PIPETTES
    left_pipette = protocol.load_instrument(
        "flex_1channel_1000", "left", tip_racks=tips
    )
    right_pipette = protocol.load_instrument(
        "flex_8channel_1000", "right", tip_racks=tips
    )
    hs_mod = protocol.load_module("heaterShakerModuleV1", "D1")
    hs_mod.open_labware_latch()
    # LOADING LABWARE
    working_reagent_reservoir = protocol.load_labware("nest_12_reservoir_15ml", "D2")
    sample_plate_slots = ["B3", "C1", "C2", "C3"]
    sample_plate = [
        protocol.load_labware("corning_96_wellplate_360ul_flat", slot)
        for slot in sample_plate_slots
    ]
    
    reagent_stock = protocol.load_labware(
        "opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical", "A1"
    )
    sample_stock_slots = ["B1", "B2"]
    sample_stock = [
        protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", slot)
        for slot in sample_stock_slots
    ]
    staging_slots = ["A4", "B4", "C4", "D4"]
    staging_racks = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_1000uL", slot)
        for slot in staging_slots
    ]


    count = 0
    samples_per_plate= 24       #ADD THIS FOR DUPLICATE AND TRIPLICATE
    num_sample_plates = math.ceil(number_samples/samples_per_plate)
    # DEFINING LIQUIDS
    bsa_stock = protocol.define_liquid(
        "BSA Stock", "BSA Stock from Pierce BCA Protein protocol ; 1.5mg/mL", "#FF6433"
    )
    bsa_rack = protocol.load_labware(
        "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap", "A2"
    )
    Reagent_A = protocol.define_liquid(
        "Reagent_A", "Reagent A for Working Reagent, will add 50 parts", "#FDF740"
    )
    Reagent_B = protocol.define_liquid(
        "Reagent_B", "Reagent B for Working Reagent, will add 1 part", "#408AFD"
    )
    water = protocol.define_liquid(
        "Diluent", "Diluent for standards, same as diluent in sample", "#D2E2FB"
    )
    sample = protocol.define_liquid(
        "sample",
        "Unknown Samples from CSV File",
        "#40FDF4",
    )

    empty_tube = protocol.define_liquid("empty", "Empty Tubes for Standards", "#D3D3D3")

    # LOADING LIQUIDS
    reagent_stock["A1"].load_liquid(water, 9000)
    bsa_rack["A1"].load_liquid(bsa_stock, 550)
    # reagent_stock["A3"].load_liquid(Reagent_A, 22000)
    # bsa_rack["D1"].load_liquid(Reagent_B, 1000)
    bsa_rack["B1"].load_liquid(empty_tube, 1)  # 1500 µg/mL
    bsa_rack["B2"].load_liquid(empty_tube, 1)  # 1000 µg/mL
    bsa_rack["B3"].load_liquid(empty_tube, 1)  # 750 µg/mL
    bsa_rack["B4"].load_liquid(empty_tube, 1)  # 500 µg/mL
    bsa_rack["B5"].load_liquid(empty_tube, 1)  # 250 µg/mL
    bsa_rack["B6"].load_liquid(empty_tube, 1)  # 125 µg/mL
    bsa_rack["C1"].load_liquid(empty_tube, 1)  # 25 µg/mL
    # bsa_rack["D6"].load_liquid(sample, 1)

    # dye_location = reagent_stock["A3"]
    dilutent_location = reagent_stock["A1"]
    sample_location = bsa_rack["D6"]
    bsa_stock_location = bsa_rack["A1"]
    
    # heatshaker.open_labware_latch()
    #Diluting Sample
    diluted_sample_offset = 6
    if protocol.params.dulute_with_walt:
        pick_up(left_pipette)
        vol_in_15_facon = get_vol_15ml_falcon(find_aspirate_height(left_pipette, dilutent_location))
        # num_transfers = math.ceil((number_samples*protocol.params.buffer_vol)/pipette_max)
        num_transfers = math.ceil((number_samples*protocol.params.buffer_vol)/(protocol.params.buffer_vol*math.floor(pipette_max/protocol.params.buffer_vol)))

        well_counter = 0
        col_num = 1#replication_mode+1     # col num for the sample_plate
        # col_num = replication_mode+1     # col num for the sample_plate
        stock_plate_num = 0
        for i in range (0, num_transfers):
            if i != 0 and i % 4 == 0:    # every 4th iteration
                remove_tip(left_pipette)
                # pick_up(left_pipette)
            if i != num_transfers-1:    # not on last iteration
                aspirate_vol = pipette_max - pipette_max%protocol.params.buffer_vol
            else:
                aspirate_vol = (number_samples*protocol.params.buffer_vol)-(pipette_max - pipette_max%protocol.params.buffer_vol)*(num_transfers-1)
            if left_pipette.has_tip == False:
                pick_up(left_pipette)
            left_pipette.blow_out(dilutent_location.top())
            left_pipette.aspirate(aspirate_vol+5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 1)
            for x in range (0, math.floor(aspirate_vol/protocol.params.buffer_vol)):
                left_pipette.dispense(protocol.params.buffer_vol, sample_stock[stock_plate_num].wells()[(well_counter%48) + 48], 0.75)
                well_counter += 1
                stock_plate_num = math.floor(well_counter/48)
            # remove_tip(left_pipette)
            vol_in_15_facon-=aspirate_vol+5
        remove_tip(left_pipette)
        for i in range (0, math.ceil(number_samples/8)):
            stock_plate_num = math.floor(i/diluted_sample_offset)
            pick_up(right_pipette)
            i = i%diluted_sample_offset
            right_pipette.aspirate(protocol.params.sample_vol, sample_stock[stock_plate_num]['A' + str(i+1)].bottom(0.1), 0.5)
            right_pipette.dispense(protocol.params.sample_vol, sample_stock[stock_plate_num]['A' + str(i+1+diluted_sample_offset)], 0.5)
            right_pipette.mix(3, protocol.params.sample_vol + protocol.params.buffer_vol-10, sample_stock[stock_plate_num]['A' + str(i+1+diluted_sample_offset)], 0.5)
            right_pipette.blow_out(sample_stock[stock_plate_num]['A' + str(i+1+diluted_sample_offset)].top())
            right_pipette.touch_tip(sample_stock[stock_plate_num]['A' + str(i+1+diluted_sample_offset)])
            right_pipette.aspirate(working_sample_vol*replication_mode+10, sample_stock[stock_plate_num]['A' + str(i+1+diluted_sample_offset)],0.5)
            for x in range (0,replication_mode):
                sample_plate_num = math.ceil((col_num*8)/(samples_per_plate*replication_mode))
                right_pipette.dispense(working_sample_vol, sample_plate[sample_plate_num-1]['A'+str( (col_num + sample_plate_num*replication_mode)%12 if (col_num + sample_plate_num*replication_mode)%12 !=0 else 12)].bottom(0.5), 0.5)
                col_num+=1
            remove_tip(right_pipette)
    def standard_loading(old, new):
        """
        old: well from sample stock
        new: row letter from sample plate
        """
        left_pipette.aspirate((working_sample_vol*replication_mode)*num_sample_plates+5, bsa_rack[old].bottom(1.5), 0.25)
        for x in range(0, num_sample_plates):
            for i in range(1, replication_mode+1):  # A1,A2,A3
                left_pipette.dispense(working_sample_vol, sample_plate[x][new + str(i)].bottom(0.1), 0.25)
        # remove_tip(left_pipette)
    # Standard Preparation  FINISH LATER
    # standard_vol_per_tube = 500#working_sample_vol*replication_mode+50
    standard_vol_per_tube = (working_sample_vol*replication_mode)*(math.ceil(number_samples/samples_per_plate))+50
    # dilutent_percentages = [0.25, 0.5, 0.625, 0.75, 0.875, 0.9375, 0.9875]
    dilutent_percentages = [1, 2/3, 1/2, 1/3, 1/6, 1/12, 1/60]
    buffer_vols =[]
    bsa_vols = []
    dilutent_pipette_vols = []
    serial_dilution_stock_tube = 3      # tubes 0 to 7
    amt_extra_in_1_over_12_tube=0
    for i in range(0, len(dilutent_percentages)):
        
        if standard_vol_per_tube*dilutent_percentages[i] < 5:
            #add extra to the serial_dilution_stock_tube duluted tube to make up for the 5ul minimum
            amt_extra_in_1_over_12_tube = (1/dilutent_percentages[serial_dilution_stock_tube])*standard_vol_per_tube*dilutent_percentages[i]
            amt_in_1_over_12_tube = standard_vol_per_tube + amt_extra_in_1_over_12_tube
            bsa_vols[serial_dilution_stock_tube] = amt_in_1_over_12_tube*dilutent_percentages[serial_dilution_stock_tube]
            buffer_vols[serial_dilution_stock_tube] =amt_in_1_over_12_tube*(1-dilutent_percentages[serial_dilution_stock_tube])
            bsa_vols.append(0)
            buffer_vols.append(standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_1_over_12_tube*(1-dilutent_percentages[i]))
        
        else:
            bsa_vols.append(standard_vol_per_tube*dilutent_percentages[i])
            buffer_vols.append(standard_vol_per_tube*(1-dilutent_percentages[i]))
    
    total_dilutent = 0
    for i in range(0, len(buffer_vols)):
        if total_dilutent + buffer_vols[i] < (pipette_max - 10):
            total_dilutent += buffer_vols[i]
            if i == len(buffer_vols) - 1:
                dilutent_pipette_vols.append(total_dilutent)
        else:
            dilutent_pipette_vols.append(total_dilutent)
            total_dilutent = buffer_vols[i]
            if i == len(buffer_vols) - 1:
                dilutent_pipette_vols.append(total_dilutent)
    # dilutent_pipette_vols.append(total_dilutent)
    
    print(buffer_vols)
    tube_spots = ["B1", "B2", "B3", "B4", "B5", "B6", "C1"]
    well_num = 0
    pick_up(left_pipette)
    if protocol.params.tip_type == 1000:        #P1000 tips
        for i in range(0, len(dilutent_pipette_vols)):
            vol_in_15_facon -= dilutent_pipette_vols[i]+10
            
            left_pipette.aspirate(
                dilutent_pipette_vols[i] + 10, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.5
            )
            amt_in_tip = dilutent_pipette_vols[i] + 10
            while amt_in_tip > 10:
                try:        # someties amt_in_tip is 10.0000000001
                    left_pipette.dispense(
                        buffer_vols[well_num],
                        bsa_rack[tube_spots[well_num]],
                        0.5,
                    )
                    amt_in_tip -= buffer_vols[well_num]
                except:
                    print(amt_in_tip)
                    break
                
                well_num += 1
    else:
        for i in range (0, len(buffer_vols)):
            vol_in_15_facon -= buffer_vols[i]
            transfer_amt = buffer_vols[i]
            for x in range (0, math.ceil(transfer_amt/pipette_max)):
                if x ==math.ceil(transfer_amt/pipette_max) - 1:
                    if buffer_vols[i]%pipette_max != 0:
                        left_pipette.aspirate(buffer_vols[i]%pipette_max, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.5)
                        left_pipette.dispense(
                                    buffer_vols[i]%pipette_max,
                                    bsa_rack[tube_spots[i]],
                                    0.5,
                                )
                    else:
                        left_pipette.aspirate(pipette_max, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.5)
                        left_pipette.dispense(
                                    pipette_max,
                                    bsa_rack[tube_spots[i]],
                                    0.5,
                                )
                else:
                    left_pipette.aspirate(pipette_max, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.5)
                    left_pipette.dispense(
                                pipette_max,
                                bsa_rack[tube_spots[i]],
                                0.5,
                            )
    remove_tip(left_pipette)
    print (bsa_vols)
    rack_order = ["B1", "B2", "B3", "B4", "B5", "B6", "C1"]
    well_order = ["A", "B", "C", "D", "E", "F", "G"]
    for i in range(0, len(bsa_vols)):
        pick_up(left_pipette)
        if bsa_vols[i] == 0:
            left_pipette.aspirate(amt_extra_in_1_over_12_tube, bsa_rack[tube_spots[serial_dilution_stock_tube]], 0.5)
            left_pipette.dispense(
                bsa_vols[i],
                bsa_rack[tube_spots[i]],
                0.5,
            )
        else:
            transfer_amt = bsa_vols[i]
            for x in range (0, math.ceil(bsa_vols[i]/pipette_max)):     
                if x == math.ceil(bsa_vols[i]/pipette_max) - 1:
                    left_pipette.aspirate(
                        transfer_amt%pipette_max,
                        bsa_stock_location,
                        0.5,
                    )
                    left_pipette.dispense(
                        transfer_amt%pipette_max,
                        bsa_rack[tube_spots[i]],
                        0.5,
                    )
                else:
                    left_pipette.aspirate(
                        pipette_max,
                        bsa_stock_location,
                        0.5,
                    )
                    left_pipette.dispense(
                        pipette_max,
                        bsa_rack[tube_spots[i]],
                        0.5,
                    )
        if standard_vol_per_tube > pipette_max:
            left_pipette.mix(3, pipette_max, bsa_rack[tube_spots[i]])
        else:
            left_pipette.mix(2, standard_vol_per_tube - 10, bsa_rack[tube_spots[i]])
        standard_loading(rack_order[i], well_order[i])
        # left_pipette.blow_out(bsa_rack[tube_spots[i]])
        remove_tip(left_pipette)

    # # Vial A
    # standard_loading("B1", "A")
    # # Vial B
    # standard_loading("B2", "B")
    # # Vial C
    # standard_loading("B3", "C")
    # # Vial D
    # standard_loading("B4", "D")
    # # Vial E
    # standard_loading("B5", "E")
    # # Vial F
    # standard_loading("B6", "F")
    # # Vial G
    # standard_loading("C1", "G")
    # Vial H: Blank
    pick_up(left_pipette)
    vol_in_15_facon-=working_sample_vol*replication_mode+5
    left_pipette.aspirate((working_sample_vol*replication_mode)*num_sample_plates+5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.25)
    for x in range (0, num_sample_plates):
        for i in range(1, replication_mode+1):  # A1,A2,A3
            left_pipette.dispense(working_sample_vol, sample_plate[x]["H" + str(i)].bottom(0.1), 0.25)
    remove_tip(left_pipette)

    # Adding Working Reagent to Plate
    num_columns = math.ceil(((math.ceil(number_samples / 8) * 8) * replication_mode) / 8) + num_sample_plates*replication_mode
    working_reagent_volume = 200
    pick_up(right_pipette)
    working_reagent_volume_amt = (num_columns*200*8) + 200#(working_reagent_volume*8*(math.ceil(number_samples/8)) + 1000)/8
    for x in range (0, num_sample_plates):
        if x == num_sample_plates-1:
            last_row = number_samples%8
            if last_row == 0:
                last_row = 12
        else:
            last_row = 12
        for i in range (1, last_row+1):
            working_reagent_volume_amt = working_reagent_volume_amt-(working_reagent_volume*8)
            # print(working_reagent_volume_amt)
            # print('A'+str(math.ceil(working_reagent_volume_amt/(10*1000))))
            right_pipette.aspirate(working_reagent_volume, working_reagent_reservoir['A'+str(math.ceil(working_reagent_volume_amt/(10*1000)))],0.75)
            right_pipette.dispense(working_reagent_volume, sample_plate[x]["A"+str(i)].top(-1), rate=0.3)
            right_pipette.blow_out(sample_plate[x]["A"+str(i)].top(-1))
            right_pipette.blow_out(sample_plate[x]["A"+str(i)].top(-1))
        #Move sample plate
        hs_mod.open_labware_latch()
        protocol.move_labware(sample_plate[x], hs_mod, use_gripper=True)
        hs_mod.close_labware_latch()
        hs_mod.set_and_wait_for_temperature(37)
        hs_mod.set_and_wait_for_shake_speed(400)

        # Shake For 30 Seconds
        if is_dry_run:
            protocol.delay(seconds=10)
        else:
            protocol.delay(minutes=0.5)
        hs_mod.deactivate_shaker()
        hs_mod.open_labware_latch()
        protocol.move_labware(sample_plate[x], new_location=protocol_api.OFF_DECK, use_gripper=False)

        
    remove_tip(right_pipette)

