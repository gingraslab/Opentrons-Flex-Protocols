# TODO: fix the volume to height function for 15ml falcon tubes so that it works for lower volumes
metadata = {
    "protocolName": "Single-plate Bradford protocol ",
    "author": "Nico To",
    "description": "Bradford for 1-24 samples in triplicate or 1-40 in duplicate.",
}
requirements = {"robotType": "Flex", "apiLevel": "2.21"}
import math
from opentrons import protocol_api
import re

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
        # return -3.33 * (volume**2) + 15.45 * volume + 9.50 - 1  # −3.33x2+15.45x+9.50
        return 0.1
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
        volume = (((height - 10.1667) / 6.41667) * 1000)
        return volume

def add_parameters(parameters):
    parameters.add_int(
        variable_name="number_samples",
        display_name="number_samples",
        description="Number of input samples.",
        default=16,
        minimum=0,
        maximum=40,
        unit="samples",
    )
    parameters.add_int(
        variable_name="diluton_amount",
        display_name="Dilution Amount",
        description="Amount to dilute by. 0 means no dilution",
        default=0,
        minimum=0,
        maximum=25,
        unit="times",

    )
    # parameters.add_int(
    #     variable_name="working_sample_vol",
    #     display_name="Working Sample Volume",
    #     description="Volume of working sample (Volume for the samples in the flat well plate)",
    #     default=5,
    #     minimum=5,
    #     maximum=50,
    #     unit="ul",
    # )
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
    add_lid = True  # protocol.params.add_lid
    working_sample_vol = 5#protocol.params.working_sample_vol
    pipette_max = 200-5

    # LOADING TIPS
    tips = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_200uL", slot)
        for slot in ["A3"]
    ]
    chute = protocol.load_waste_chute()

    def remove_tip(pipette):
        if is_dry_run:
            pipette.return_tip()
        else:
            pipette.drop_tip(chute)

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

    # LOADING LABWARE
    working_reagent_reservoir = protocol.load_labware("nest_12_reservoir_15ml", "D2")
    heatshaker = protocol.load_module("heaterShakerModuleV1", "D1")
    working_plate = protocol.load_labware("corning_96_wellplate_360ul_flat", "C2")
    reagent_stock = protocol.load_labware(
        "opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical", "A1"
    )
    if add_lid:
        lid = protocol.load_labware(
            "opentrons_tough_pcr_auto_sealing_lid", location="C1"
        )
    sample_stock = protocol.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", "B2"
    )
    staging_slots = ["A4", "B4", "C4"]
    staging_racks = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_1000uL", slot)
        for slot in staging_slots
    ]

    # REPLENISHING TIPS

    count = 0
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
    dilutent = protocol.define_liquid(
        "Diluent", "Diluent for standards, same as diluent in sample", "#D2E2FB"
    )
    sample = protocol.define_liquid(
        "sample",
        "sample",
        "#40FDF4",
    )

    empty_tube = protocol.define_liquid("empty", "Empty Tubes for Standards", "#D3D3D3")

    # LOADING LIQUIDS
    reagent_stock["A1"].load_liquid(dilutent, 9000)
    reagent_stock["A2"].load_liquid(Reagent_A, 9000)
    bsa_rack["B1"].load_liquid(bsa_stock, 550)
    # reagent_stock["A3"].load_liquid(Reagent_A, 22000)
    # bsa_rack["D1"].load_liquid(Reagent_B, 1000)
    # bsa_rack["B1"].load_liquid(empty_tube, 1)  # 1500 µg/mL
    bsa_rack["B2"].load_liquid(empty_tube, 1)  # 1000 µg/mL
    bsa_rack["B3"].load_liquid(empty_tube, 1)  # 750 µg/mL
    bsa_rack["B4"].load_liquid(empty_tube, 1)  # 500 µg/mL
    bsa_rack["B5"].load_liquid(empty_tube, 1)  # 250 µg/mL
    # bsa_rack["B6"].load_liquid(empty_tube, 1)  # 125 µg/mL
    # bsa_rack["C1"].load_liquid(empty_tube, 1)  # 25 µg/mL
    # bsa_rack["D6"].load_liquid(sample, 1)

    # dye_location = reagent_stock["A3"]
    dilutent_location = reagent_stock["A1"]
    reagent_a_location = reagent_stock["A2"]
    # sample_location = bsa_rack["D6"]
    bsa_stock_location = bsa_rack["B1"]
    
    #Variables for creating standards
    standard_vol_per_tube = 200#(working_sample_vol*replication_mode+50)*2
    concentrations = [1.5, 1, 0.75, 0.5, 0.25]
    tube_spots = ["B1", "B2", "B3", "B4", "B5", "B6", "C1"]
    well_order = ["A", "B", "C", "D", "E", "F", "G", "H"]
    
    heatshaker.open_labware_latch()
    
    # Adding 25 ul Reagent A to Plate
    number_occupied_wells = number_samples * replication_mode + replication_mode*(len(concentrations)+1)  # number of wells occupied by samples and standards
    amt_reagent_a = 25
    num_transfers = math.ceil((number_occupied_wells*amt_reagent_a)/(amt_reagent_a*(math.floor(pipette_max/amt_reagent_a))))
    well_counter = 0
    left_pipette.pick_up_tip()
    vol_in_15_falcon_reagent_a = get_vol_15ml_falcon(find_aspirate_height(left_pipette, reagent_a_location))




    # Wells with reagent A
    regA_occupied_wells = []
    #Standard Wells
    for i in range (0, len(concentrations)+1):
        for x in range (1, replication_mode+1):  # A1, A2, A3
            regA_occupied_wells.append(well_order[i] + str(x))  # A1, A2, A3, B1, B2, B3, C1
    current_column = replication_mode+ 1
    for i in range (0, math.ceil(number_samples/8)):
        if i != math.ceil(number_samples/8)-1:  # not on last iteration
            for letter in well_order:
                for x in range (1, replication_mode+1):
                    regA_occupied_wells.append(letter + str(current_column+x-1))
            current_column += replication_mode
            print("1a")
        else:
            print("2")
            remainder = number_samples % 8 if number_samples % 8 != 0 else 8
            for letter in well_order[0:remainder]:
                for x in range (1, replication_mode+1):
                    regA_occupied_wells.append(letter + str(current_column+x-1))
            current_column += replication_mode
    regA_occupied_wells = sorted(regA_occupied_wells)
    regA_occupied_wells = sorted(regA_occupied_wells, key=lambda x: (x[0], int(re.findall(r'\d+', x)[0])))

    print(regA_occupied_wells)
    
    for i in range (0, num_transfers):      # FIX THIS
        if i != num_transfers-1:    # not on last iteration
            aspirate_vol = pipette_max - pipette_max%amt_reagent_a
        else:
            aspirate_vol = (number_occupied_wells*amt_reagent_a)-(pipette_max - pipette_max%amt_reagent_a)*(num_transfers-1)
        # print(aspirate_vol)
        if left_pipette.has_tip == False:
            left_pipette.pick_up_tip()
        left_pipette.blow_out(reagent_a_location.top())
        try:
            left_pipette.aspirate(aspirate_vol+5, reagent_a_location.bottom(get_height_15ml_falcon(vol_in_15_falcon_reagent_a)), 0.5)
        except:
            left_pipette.aspirate(aspirate_vol+5, reagent_a_location.bottom(1), 0.5)

        # left_pipette.aspirate(aspirate_vol+5, reagent_a_location.bottom(1), 0.5)
        for x in range (0, math.floor(aspirate_vol/amt_reagent_a)):
            left_pipette.dispense(amt_reagent_a, working_plate[regA_occupied_wells[well_counter]].bottom(0.2), 0.1)
            well_counter += 1
        # remove_tip(left_pipette)
        vol_in_15_falcon_reagent_a-=aspirate_vol+5
        remove_tip(left_pipette)
    # remove_tip(left_pipette)

    
    
    #Diluting Sample
    if protocol.params.diluton_amount == 0:
        dilute_with_walt = False
    else:
        dilute_with_walt = True
    
    if dilute_with_walt:
        sample_vol = max((working_sample_vol*3+5)/protocol.params.diluton_amount, 5)
        
        for i in range (0, number_samples):
            sample_stock.wells()[i].load_liquid(sample, sample_vol)

        buffer_vol = sample_vol*protocol.params.diluton_amount - sample_vol
        diluted_sample_offset = 6
        left_pipette.pick_up_tip()
        vol_in_15_falcon_dilutent =  get_vol_15ml_falcon(find_aspirate_height(left_pipette, dilutent_location))
        num_transfers = math.ceil((number_samples*buffer_vol)/pipette_max)
        well_counter = 0
        col_num = replication_mode+1     # col num for the working_plate
        for i in range (0, num_transfers):
            if left_pipette.has_tip == False:
                left_pipette.pick_up_tip()
            
            if i != num_transfers-1:    # not on last iteration
                aspirate_vol = min(pipette_max - pipette_max%buffer_vol, buffer_vol)
            else:
                aspirate_vol = min((number_samples*buffer_vol)-(pipette_max - pipette_max%buffer_vol)*(num_transfers-1), buffer_vol)
            if left_pipette.has_tip == False:
                left_pipette.pick_up_tip()
            left_pipette.blow_out(dilutent_location.top())
            try:
                left_pipette.aspirate(aspirate_vol+5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_falcon_dilutent)), 1)
            except:
                left_pipette.aspirate(aspirate_vol+5, dilutent_location.bottom(1), 1)

            for x in range (0, math.floor(aspirate_vol/buffer_vol)):
                left_pipette.dispense(buffer_vol, sample_stock.wells()[well_counter + 48], 0.75)
                well_counter += 1
            remove_tip(left_pipette)
            vol_in_15_falcon_dilutent-=aspirate_vol+5
            # if i %3 == 0 and i != 0:
            #     remove_tip(left_pipette)
        if left_pipette.has_tip:
            remove_tip(left_pipette)
        for i in range (0, math.ceil(number_samples/8)):
            right_pipette.pick_up_tip()
            print(sample_vol)
            right_pipette.aspirate(sample_vol, sample_stock['A' + str(i+1)].bottom(0.1), 0.1)
            right_pipette.dispense(sample_vol, sample_stock['A' + str(i+1+diluted_sample_offset)], 0.1)
            right_pipette.mix(3, sample_vol + buffer_vol-5, sample_stock['A' + str(i+1+diluted_sample_offset)], 0.1)
            right_pipette.blow_out(sample_stock['A' + str(i+1+diluted_sample_offset)].top())
            right_pipette.touch_tip(sample_stock['A' + str(i+1+diluted_sample_offset)])
            right_pipette.aspirate(working_sample_vol*replication_mode+10, sample_stock['A' + str(i+1+diluted_sample_offset)],0.1)
            for x in range (0,replication_mode):
                right_pipette.dispense(working_sample_vol, working_plate['A' + str(col_num)].bottom(0.2), 0.1)
                # right_pipette.blow_out(working_plate['A' + str(col_num)].top())
                col_num+=1
            remove_tip(right_pipette)
    else:
        for i in range (0, number_samples):
            sample_stock.wells()[i].load_liquid(sample, working_sample_vol*3*10)

        col_num = replication_mode+1
        for i in range (0, math.ceil(number_samples/8)):
            right_pipette.pick_up_tip()
            right_pipette.aspirate(working_sample_vol*3+5, sample_stock['A' + str(i+1)],0.1)
            for x in range (0,replication_mode):
                right_pipette.dispense(working_sample_vol, working_plate['A' + str(col_num)].bottom(0.2), 0.1)
                # right_pipette.blow_out(working_plate['A' + str(col_num)].top())
                col_num+=1
            remove_tip(right_pipette)
    def standard_loading(old, new):
        """
        old: well from sample stock
        new: row letter from sample plate
        """

        if left_pipette.has_tip == False:
            left_pipette.pick_up_tip()
        left_pipette.blow_out(bsa_rack[old].top())
        remove_tip(left_pipette)
        left_pipette.pick_up_tip()
        
        for i in range(1, replication_mode+1):  # A1,A2,A3
            left_pipette.aspirate(working_sample_vol, bsa_rack[old].bottom(1.5), 0.1)
            left_pipette.dispense(working_sample_vol, working_plate[new + str(i)].bottom(0.2), 0.1)
            left_pipette.blow_out(working_plate[new + str(i)].top(-5))
            left_pipette.touch_tip(working_plate[new + str(i)])
        # remove_tip(left_pipette)

        
        # left_pipette.aspirate(working_sample_vol*3+5, bsa_rack[old].bottom(1.5), 0.1)
        # for i in range(1, replication_mode+1):  # A1,A2,A3
        #     left_pipette.dispense(working_sample_vol, working_plate[new + str(i)].bottom(0.2), 0.2)
        #     # left_pipette.blow_out(working_plate[new + str(i)].top())
        # remove_tip(left_pipette)


    #Loading buffer for standard preparation
    buffer_amts = []
    for i in range(0, len(concentrations)):
        if concentrations[i] == 1.5:
            buffer_amts.append(0)
            continue
        else:       # serial dilution with previous tube
            print(standard_vol_per_tube)
            print(concentrations[i])
            amt_extra_in_tube = (concentrations[i]*standard_vol_per_tube)/concentrations[i-1]#(1/dilutent_percentages[serial_dilution_stock_tube])*standard_vol_per_tube*dilutent_percentages[i]
            # print(amt_extra_in_tube)
            # buffer_amt = standard_vol_per_tube-amt_extra_in_1_over_12_tube#standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_tube*(1-dilutent_percentages[i])
            buffer_amt = standard_vol_per_tube - amt_extra_in_tube#standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_tube*(1-dilutent_percentages[i])
            buffer_amts.append(buffer_amt)
    num_transfers = math.ceil(sum(buffer_amts)/pipette_max)
    transfers = []
    temp_amt = 0
    for i in range (0, len(buffer_amts)):
        if temp_amt + buffer_amts[i] > pipette_max:
            transfers.append(temp_amt)
            temp_amt = buffer_amts[i]
            if i == (len(buffer_amts)-1):
                transfers.append(temp_amt)

        else:
            temp_amt += buffer_amts[i]
            if i == (len(buffer_amts)-1):
                transfers.append(temp_amt)
        # print(temp_amt)
    tube_tracker = 0
    # print(transfers)
    left_pipette.pick_up_tip()
    vol_in_15_falcon_dilutent= get_vol_15ml_falcon(find_aspirate_height(left_pipette, dilutent_location))
    for i in range(0, num_transfers):
        # left_pipette.pick_up_tip()
        vol_in_15_falcon_dilutent -= buffer_amt
        left_pipette.blow_out(dilutent_location.top())
        try:
            left_pipette.aspirate(transfers[i] +5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_falcon_dilutent)), 0.1)
        except:
            left_pipette.aspirate(transfers[i] +5, dilutent_location.bottom(1), 0.1)
        while transfers[i] > 0:
            if buffer_amts[tube_tracker] != 0:
                left_pipette.dispense(buffer_amts[tube_tracker], bsa_rack[tube_spots[tube_tracker]], rate=0.1)
                # print(buffer_amts[tube_tracker])
            transfers[i] -= buffer_amts[tube_tracker]
            tube_tracker += 1
    remove_tip(left_pipette)
    print(buffer_amts)
    
    # Standard Preparation  
    for i in range(0, len(concentrations)):
        if concentrations[i] == 1.5:
            #buffer
            # buffer_vol = (standard_vol_per_tube)*(1-dilutent_percentages[i])
            if left_pipette.has_tip == False:
                left_pipette.pick_up_tip()
            # remove_tip(left_pipette)
            standard_loading(tube_spots[i], well_order[i])   
            # remove_tip(left_pipette)     
            continue

        else:       # serial dilution with previous tube
            print(standard_vol_per_tube)
            print(concentrations[i])
            amt_extra_in_tube = (concentrations[i]*standard_vol_per_tube)/concentrations[i-1]#(1/dilutent_percentages[serial_dilution_stock_tube])*standard_vol_per_tube*dilutent_percentages[i]
            # print(amt_extra_in_tube)
            # buffer_amt = standard_vol_per_tube-amt_extra_in_1_over_12_tube#standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_tube*(1-dilutent_percentages[i])
            buffer_amt = standard_vol_per_tube - amt_extra_in_tube#standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_tube*(1-dilutent_percentages[i])
            
            # #buffer
            if left_pipette.has_tip == False:
                left_pipette.pick_up_tip()
            # vol_in_15_falcon_dilutent -= buffer_amt
            # print("current : " + str(concentrations[i]))
            # print("buffer: " + str(buffer_amt))
            # left_pipette.aspirate(buffer_amt, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_falcon_dilutent)), 0.1)
            # left_pipette.dispense(buffer_amt, bsa_rack[tube_spots[i]], rate=0.2)
            # left_pipette.blow_out(bsa_rack[tube_spots[i]].top())
            # remove_tip(left_pipette)
            #bsa
            # left_pipette.pick_up_tip()
            print("Stock: " + str(amt_extra_in_tube))
            left_pipette.aspirate(amt_extra_in_tube, bsa_rack[tube_spots[i-1]], 0.1)
            left_pipette.dispense(amt_extra_in_tube, bsa_rack[tube_spots[i]], 0.1)
            left_pipette.mix(3, standard_vol_per_tube-5, bsa_rack[tube_spots[i]], 0.3)
            left_pipette.blow_out(bsa_rack[tube_spots[i]].top(1))
            standard_loading(tube_spots[i], well_order[i])
            # remove_tip(left_pipette)
    
   
    # try:
    #     vol_in_15_falcon_dilutent
    # except NameError:
    #     vol_in_15_falcon_dilutent = get_vol_15ml_falcon(find_aspirate_height(left_pipette, dilutent_location))

    vol_in_15_falcon_dilutent-=working_sample_vol*replication_mode
    remove_tip(left_pipette)
    left_pipette.pick_up_tip()
    try:
        left_pipette.aspirate(working_sample_vol*replication_mode+5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_falcon_dilutent)), 0.25)
    except:
        left_pipette.aspirate(working_sample_vol*replication_mode+5, dilutent_location.bottom(1), 0.25)
    for i in range(1, replication_mode+1):  # A1,A2,A3
        current_letter = well_order[len(concentrations)]
        left_pipette.dispense(working_sample_vol, working_plate[current_letter + str(i)].bottom(0.1), 0.1)
    remove_tip(left_pipette)

    # Adding Working Reagent (Reagent B) to Plate
    num_columns = math.ceil(((math.ceil(number_samples / 8) * 8) * replication_mode) / 8) + replication_mode
    working_reagent_volume = 200
    right_pipette.pick_up_tip()
    working_reagent_volume_amt = num_columns*200*8 +1000#(working_reagent_volume*8*(math.ceil(number_samples/8)) + 1000)/8
    dye = protocol.define_liquid(
        "Dye", "Dye", "#A840FD"
    )

    for i in range (0, math.ceil(working_reagent_volume_amt/(10.5*1000))):
        if i ==math.ceil(working_reagent_volume_amt/(10.5*1000))-1:
            working_reagent_reservoir.wells()[i].load_liquid(dye, working_reagent_volume_amt-10500*i + 1500)
        else:
            working_reagent_reservoir.wells()[i].load_liquid(dye, 12000)

    
    for i in range (0, num_columns):
        working_reagent_volume_amt = working_reagent_volume_amt-(working_reagent_volume*8)
        # protocol.comment(str(working_reagent_volume_amt))
        # protocol.comment(str(math.ceil(working_reagent_volume_amt/(10.5*1000))))
        right_pipette.aspirate(working_reagent_volume, working_reagent_reservoir['A'+str(math.ceil(working_reagent_volume_amt/(10.5*1000)))],0.5)
        right_pipette.dispense(working_reagent_volume, working_plate["A"+str(i+1)].top(-1), rate=0.3)
        right_pipette.blow_out(working_plate["A"+str(i+1)].top(-1))
        right_pipette.blow_out(working_plate["A"+str(i+1)].top(-1))
    remove_tip(right_pipette)

    # Prep HeaterShaker
    heatshaker.open_labware_latch()
    # move labware with lid onto hs_mod
    if add_lid:
        protocol.deck.__delitem__("C2")
        new_working_plate = protocol.load_labware(
            "opentrons_96_wellplate_200ul_pcr_full_skirt", "C2"
        )
        protocol.move_labware(
            new_working_plate, new_location=heatshaker, use_gripper=True
        )
        heatshaker.close_labware_latch()
        heatshaker.open_labware_latch()
        new_working_plate.set_offset(x=0.00, y=0.00, z=30)
        protocol.move_labware(
            labware=lid, new_location=new_working_plate, use_gripper=True
        )
    else:
        protocol.move_labware(
            labware=working_plate, new_location=heatshaker, use_gripper=True
        )
    # protocol.pause("Place lid on well plate")
    heatshaker.close_labware_latch()
    # heatshaker.set_and_wait_for_temperature(37)
    heatshaker.set_and_wait_for_shake_speed(400)

    # Shake For 30 Seconds
    if is_dry_run:
        protocol.delay(seconds=10)
    else:
        protocol.delay(minutes=0.5)
    heatshaker.deactivate_shaker()
    heatshaker.open_labware_latch()

    protocol.comment("\n---------------15 Minute Incubation----------------\n\n")
    if is_dry_run:
        protocol.delay(seconds=10)
    else:
        protocol.delay(minutes=15)  # SEND EMAIL AT 10 MINUTES

    # Deactivating Heatshaker
    heatshaker.deactivate_heater()
    heatshaker.open_labware_latch()
    if add_lid:
        protocol.move_labware(lid, "C1", use_gripper=True)
        protocol.move_labware(new_working_plate, "C2", use_gripper=True)
        heatshaker.close_labware_latch()
    else:
        protocol.move_labware(working_plate, "C2", use_gripper=True)
        heatshaker.close_labware_latch()
    # left_pipette.pick_up_tip()
    # left_pipette.pick_up_tip()